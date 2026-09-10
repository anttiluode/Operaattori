from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "audits"
if str(AUDITS) not in sys.path:
    sys.path.insert(0, str(AUDITS))

from cross_cell_gradient_decomposition import site_directions
from cross_cell_operator import (
    FCI_COMMIT,
    SITE_X,
    choose_branches,
    git_head,
    instantiate_matched_passive,
    load_panel,
    section_name,
    setup_neuron,
)
from direct_cable_graph import build_compartment_graph, node_for_section_x
from gradient_operating_point import conductance_program
from gradient_transport_feedback import simulate_decomposition
from real_metric_tangent import assemble_metric_graph, compile_sparse, simulate

FINITE_SIGN_EPS_MV = 1.0e-5
TANGENT_SIGN_EPS = 1.0e-12
ALLOWED_DRIVES = (0.5, 1.0, 2.0)
ALLOWED_Q = (1.05, 1.10)


def rms(x: np.ndarray) -> float:
    a = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(a * a)))


def relerr(actual: np.ndarray, predicted: np.ndarray) -> float:
    return rms(np.asarray(predicted) - np.asarray(actual)) / (rms(actual) + 1e-30)


def tangent_sign(value: float) -> int:
    if value > TANGENT_SIGN_EPS:
        return 1
    if value < -TANGENT_SIGN_EPS:
        return -1
    return 0


def finite_sign(value: float) -> int:
    if value > FINITE_SIGN_EPS_MV:
        return 1
    if value < -FINITE_SIGN_EPS_MV:
        return -1
    return 0


def replay_frozen_current(
    G,
    C: np.ndarray,
    site_nodes: np.ndarray,
    soma_node: int,
    current_uA: np.ndarray,
) -> np.ndarray:
    """Replay a fixed site-current waveform through a compiled passive graph."""
    comp = compile_sparse(G, C, site_nodes, soma_node)
    d = comp["d"]
    lu = comp["lu"]
    X = comp["X"]
    state = np.zeros(len(C), dtype=float)
    soma = np.zeros(current_uA.shape[1], dtype=float)
    for ti in range(current_uA.shape[1]):
        passive = lu.solve(d * state)
        state = passive + X @ current_uA[:, ti]
        soma[ti] = state[int(soma_node)]
    return soma


def summarize_records(records: list[dict]) -> dict:
    full_eval = [r for r in records if r["finite_full_sign"] != 0]
    transport_eval = [r for r in records if r["finite_transport_sign"] != 0]
    predicted_reversal_eval = [
        r
        for r in records
        if r["predicted_reversal"]
        and r["finite_full_sign"] != 0
        and r["finite_transport_sign"] != 0
    ]
    actual_reversals = [r for r in records if r["actual_finite_reversal"]]
    true_positive = [
        r for r in records if r["predicted_reversal"] and r["actual_finite_reversal"]
    ]

    def accuracy(rows: list[dict], predicted_key: str, actual_key: str) -> float | None:
        if not rows:
            return None
        return float(np.mean([r[predicted_key] == r[actual_key] for r in rows]))

    precision = (
        float(len(true_positive) / len(predicted_reversal_eval))
        if predicted_reversal_eval
        else None
    )
    recall = (
        float(len(true_positive) / len(actual_reversals))
        if actual_reversals
        else None
    )

    return {
        "records": len(records),
        "full_sign_evaluable": len(full_eval),
        "transport_sign_evaluable": len(transport_eval),
        "full_tangent_finite_sign_accuracy": accuracy(
            full_eval, "tangent_full_sign", "finite_full_sign"
        ),
        "transport_tangent_finite_sign_accuracy": accuracy(
            transport_eval, "tangent_transport_sign", "finite_transport_sign"
        ),
        "predicted_reversals": int(sum(r["predicted_reversal"] for r in records)),
        "predicted_reversals_finite_evaluable": len(predicted_reversal_eval),
        "actual_finite_reversals": len(actual_reversals),
        "reversal_true_positive": len(true_positive),
        "predicted_reversal_precision": precision,
        "predicted_reversal_recall": recall,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--drive-scale", type=float, required=True)
    ap.add_argument("--metric-scale", type=float, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    drive = float(args.drive_scale)
    q = float(args.metric_scale)
    if drive not in ALLOWED_DRIVES:
        raise ValueError(f"drive must be one of {ALLOWED_DRIVES}")
    if q not in ALLOWED_Q:
        raise ValueError(f"metric scale must be one of {ALLOWED_Q}")
    log_step = float(math.log(q))

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source not pinned")

    setup_neuron(fci_root)
    panel = load_panel(fci_root)
    ga, gn = conductance_program(drive)
    records: list[dict] = []
    operating_points: list[dict] = []

    for ci, model in enumerate(panel):
        cell = instantiate_matched_passive(fci_root, model)
        graph = build_compartment_graph(cell)
        branch = choose_branches(cell)[0]
        site_nodes = np.asarray(
            [node_for_section_x(graph, branch, float(x)) for x in SITE_X], dtype=int
        )
        soma_node = int(node_for_section_x(graph, cell.soma[0], 0.5))
        directions, meta = site_directions(graph, branch)

        base = simulate_decomposition(
            graph["G_uS"],
            graph["C_nF"],
            site_nodes,
            soma_node,
            ga,
            gn,
            directions,
        )
        if not base["all_converged"]:
            raise RuntimeError(
                f"base nonlinear solve failed for {model['morphology_identifier']}"
            )
        peak_index = int(np.argmax(base["soma_mV"]))
        base_fixed = float(base["soma_mV"][peak_index])
        full_tangent = np.asarray(base["full_tangent"][:, peak_index], dtype=float)
        transport_tangent = np.asarray(
            base["frozen_current_tangent"][:, peak_index], dtype=float
        )
        base_current = np.asarray(base["current_uA"], dtype=float)

        replay_base = replay_frozen_current(
            graph["G_uS"],
            graph["C_nF"],
            site_nodes,
            soma_node,
            base_current,
        )
        replay_identity_error = relerr(base["soma_mV"], replay_base)
        operating_points.append(
            {
                "cell_order": int(model["order"]),
                "morphology_identifier": model["morphology_identifier"],
                "species": model["species"],
                "layer": model["layer"],
                "compartments": int(len(graph["nodes"])),
                "branch": section_name(branch),
                "drive_scale": drive,
                "base_peak_mV": base_fixed,
                "base_peak_index": peak_index,
                "base_peak_time_ms": float((peak_index + 1) * 0.025),
                "base_replay_relative_rms_error": replay_identity_error,
            }
        )

        for di, (mi, gf, gt) in enumerate(zip(meta, full_tangent, transport_tangent)):
            Gq, Cq = assemble_metric_graph(
                graph,
                target_node=int(mi["node"]),
                kind=str(mi["kind"]),
                log_scale=log_step,
            )
            full_out = simulate(Gq, Cq, site_nodes, soma_node, ga, gn)
            replay_out = replay_frozen_current(
                Gq, Cq, site_nodes, soma_node, base_current
            )
            delta_full = float(full_out["soma_mV"][peak_index] - base_fixed)
            delta_transport = float(replay_out[peak_index] - base_fixed)
            pred_full = float(gf * log_step)
            pred_transport = float(gt * log_step)
            s_gf = tangent_sign(float(gf))
            s_gt = tangent_sign(float(gt))
            s_ff = finite_sign(delta_full)
            s_ft = finite_sign(delta_transport)
            records.append(
                {
                    "cell_order": int(model["order"]),
                    "morphology_identifier": model["morphology_identifier"],
                    "species": model["species"],
                    "layer": model["layer"],
                    "branch": section_name(branch),
                    "drive_scale": drive,
                    "metric_scale": q,
                    "log_metric_step": log_step,
                    "direction_index": int(di),
                    **mi,
                    "base_peak_mV": base_fixed,
                    "full_tangent_mV_per_logscale": float(gf),
                    "transport_tangent_mV_per_logscale": float(gt),
                    "feedback_tangent_mV_per_logscale": float(gf - gt),
                    "predicted_full_change_mV": pred_full,
                    "predicted_transport_change_mV": pred_transport,
                    "finite_full_change_mV": delta_full,
                    "finite_transport_change_mV": delta_transport,
                    "tangent_full_sign": s_gf,
                    "tangent_transport_sign": s_gt,
                    "finite_full_sign": s_ff,
                    "finite_transport_sign": s_ft,
                    "predicted_reversal": bool(
                        s_gf != 0 and s_gt != 0 and s_gf != s_gt
                    ),
                    "actual_finite_reversal": bool(
                        s_ff != 0 and s_ft != 0 and s_ff != s_ft
                    ),
                    "full_all_converged": bool(full_out["all_converged"]),
                    "full_abs_prediction_error_mV": float(abs(pred_full - delta_full)),
                    "transport_abs_prediction_error_mV": float(
                        abs(pred_transport - delta_transport)
                    ),
                    "full_relative_prediction_error": (
                        float(abs(pred_full - delta_full) / abs(delta_full))
                        if abs(delta_full) >= FINITE_SIGN_EPS_MV
                        else None
                    ),
                    "transport_relative_prediction_error": (
                        float(abs(pred_transport - delta_transport) / abs(delta_transport))
                        if abs(delta_transport) >= FINITE_SIGN_EPS_MV
                        else None
                    ),
                }
            )

        print(
            f"[{ci+1:02d}/24] drive={drive:g} q={q:g} "
            f"{model['species']:5s} {model['morphology_identifier']:>12s} "
            f"replay={replay_identity_error:.3e}",
            flush=True,
        )

    result = {
        "object": "finite prediction of transport versus NMDA-feedback morphology derivative split",
        "protocol": "FINITE_DECOMPOSITION_PROTOCOL.md",
        "fci_commit": FCI_COMMIT,
        "drive_scale": drive,
        "metric_scale": q,
        "log_metric_step": log_step,
        "finite_sign_epsilon_mV": FINITE_SIGN_EPS_MV,
        "operating_points": operating_points,
        "records": records,
        "summary": summarize_records(records),
        "max_base_replay_relative_rms_error": float(
            max(x["base_replay_relative_rms_error"] for x in operating_points)
        ),
        "all_finite_full_solves_converged": bool(
            all(r["full_all_converged"] for r in records)
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
