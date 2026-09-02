from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "audits"
if str(AUDITS) not in sys.path:
    sys.path.insert(0, str(AUDITS))

from causal_nonlinear_graph import (
    DT_MS,
    REST_MV,
    current_derivative,
    solve_site_step,
    synapse_law,
)
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
from direct_cable_graph import (
    build_compartment_graph,
    node_for_section_x,
)
from gradient_operating_point import (
    DRIVE_SCALES,
    branch_directions,
    conductance_program,
)
from real_metric_tangent import (
    compile_sparse,
    find_cell1125,
    simulate_with_metric_tangents,
)


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def relerr(actual: np.ndarray, pred: np.ndarray) -> float:
    return rms(np.asarray(pred) - np.asarray(actual)) / (
        rms(actual) + 1e-30
    )


def sign(v: float, eps: float = 1e-12) -> int:
    if v > eps:
        return 1
    if v < -eps:
        return -1
    return 0


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(
        np.dot(a, b)
        / ((np.linalg.norm(a) * np.linalg.norm(b)) + 1e-30)
    )


def simulate_decomposition(
    G,
    C,
    site_nodes: np.ndarray,
    soma_node: int,
    ga: np.ndarray,
    gn: np.ndarray,
    directions,
) -> dict:
    """Split dV/dtheta into frozen-current transport + NMDA feedback.

    The base nonlinear trajectory is unchanged.

    frozen-current transport tangent:
        differentiate the compiled passive recurrence and X(theta) while
        replaying the *base synaptic current waveform* as an external input.

    full tangent:
        additionally differentiate J(V), including the same local implicit
        Jacobian used by the causal NMDA solve.

    feedback contribution is defined exactly as full - frozen.
    """
    comp = compile_sparse(G, C, site_nodes, soma_node)
    d = comp["d"]
    lu = comp["lu"]
    X = comp["X"]
    R = comp["R"]

    n = len(C)
    k = len(directions)
    dd = []
    dA = []
    dX = []
    dR = []
    for dG, dC in directions:
        ddi = np.asarray(dC, dtype=float) / DT_MS
        dAi = (
            dG
            + sparse.diags(
                ddi,
                offsets=0,
                shape=(n, n),
                format="csr",
            )
        ).tocsr()
        dXi = lu.solve(-(dAi @ X))
        dd.append(ddi)
        dA.append(dAi)
        dX.append(dXi)
        dR.append(dXi[site_nodes, :])

    state = np.zeros(n, dtype=float)
    dfull = np.zeros((n, k), dtype=float)
    dfrozen = np.zeros((n, k), dtype=float)
    previous_z = np.zeros(len(site_nodes), dtype=float)

    ntime = ga.shape[1]
    soma = np.zeros(ntime, dtype=float)
    full_trace = np.zeros((k, ntime), dtype=float)
    frozen_trace = np.zeros((k, ntime), dtype=float)
    current_trace = np.zeros((len(site_nodes), ntime), dtype=float)
    local_trace = np.zeros((len(site_nodes), ntime), dtype=float)

    all_converged = True
    max_newton = 0

    for ti in range(ntime):
        old_state = state
        old_full = dfull
        old_frozen = dfrozen

        passive = lu.solve(d * old_state)

        common = np.empty((n, k), dtype=float)
        for pi in range(k):
            common[:, pi] = (
                dd[pi] * old_state
                - dA[pi] @ passive
            )

        rhs_full = common + d[:, None] * old_full
        rhs_frozen = common + d[:, None] * old_frozen
        dpassive_full = lu.solve(rhs_full)
        dpassive_frozen = lu.solve(rhs_frozen)

        solved = solve_site_step(
            passive[site_nodes],
            R,
            ga[:, ti],
            gn[:, ti],
            previous_z,
        )
        all_converged = (
            all_converged and bool(solved["converged"])
        )
        max_newton = max(
            max_newton, int(solved["iterations"])
        )

        z = np.asarray(solved["z"], dtype=float)
        absolute_v = REST_MV + z
        current = synapse_law(
            absolute_v,
            ga[:, ti],
            gn[:, ti],
        )
        jv = current_derivative(
            absolute_v,
            ga[:, ti],
            gn[:, ti],
        )

        K = np.eye(len(site_nodes)) - R @ np.diag(jv)
        dz_rhs = np.empty((len(site_nodes), k), dtype=float)
        for pi in range(k):
            dz_rhs[:, pi] = (
                dpassive_full[site_nodes, pi]
                + dR[pi] @ current
            )
        dz = np.linalg.solve(K, dz_rhs)
        dcurrent = jv[:, None] * dz

        state = passive + X @ current
        dfull = np.empty_like(old_full)
        dfrozen = np.empty_like(old_frozen)
        for pi in range(k):
            geometric_injection = dX[pi] @ current
            dfull[:, pi] = (
                dpassive_full[:, pi]
                + geometric_injection
                + X @ dcurrent[:, pi]
            )
            dfrozen[:, pi] = (
                dpassive_frozen[:, pi]
                + geometric_injection
            )

        soma[ti] = state[soma_node]
        full_trace[:, ti] = dfull[soma_node, :]
        frozen_trace[:, ti] = dfrozen[soma_node, :]
        current_trace[:, ti] = current
        local_trace[:, ti] = state[site_nodes]
        previous_z = state[site_nodes].copy()

    return {
        "soma_mV": soma,
        "full_tangent": full_trace,
        "frozen_current_tangent": frozen_trace,
        "feedback_tangent": full_trace - frozen_trace,
        "current_uA": current_trace,
        "local_depol_mV": local_trace,
        "all_converged": bool(all_converged),
        "max_newton_iterations": int(max_newton),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/"
            "gradient_transport_feedback.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source not pinned")

    setup_neuron(fci_root)
    model = find_cell1125(load_panel(fci_root))
    cell = instantiate_matched_passive(fci_root, model)
    graph = build_compartment_graph(cell)
    branch = choose_branches(cell)[0]

    site_nodes = np.asarray(
        [
            node_for_section_x(graph, branch, float(x))
            for x in SITE_X
        ],
        dtype=int,
    )
    soma_node = node_for_section_x(graph, cell.soma[0], 0.5)
    directions, meta = branch_directions(graph, branch)

    records = []
    reference_error = None

    for scale in DRIVE_SCALES:
        ga, gn = conductance_program(scale)
        result = simulate_decomposition(
            graph["G_uS"],
            graph["C_nF"],
            site_nodes,
            soma_node,
            ga,
            gn,
            directions,
        )
        peak_index = int(np.argmax(result["soma_mV"]))
        full = result["full_tangent"][:, peak_index]
        frozen = result[
            "frozen_current_tangent"
        ][:, peak_index]
        feedback = result["feedback_tangent"][:, peak_index]

        if float(scale) == 1.0:
            reference = simulate_with_metric_tangents(
                graph["G_uS"],
                graph["C_nF"],
                site_nodes,
                soma_node,
                ga,
                gn,
                directions,
            )
            reference_error = relerr(
                reference["soma_tangent_mV_per_logscale"],
                result["full_tangent"],
            )

        rows = []
        for mi, gf, gt, gnlf in zip(
            meta, full, frozen, feedback
        ):
            rows.append(
                {
                    **mi,
                    "full_dpeak_mV_per_logscale": float(gf),
                    "frozen_transport_dpeak_mV_per_logscale": float(gt),
                    "nmda_feedback_dpeak_mV_per_logscale": float(gnlf),
                    "full_sign": sign(gf),
                    "frozen_sign": sign(gt),
                    "feedback_sign": sign(gnlf),
                    "feedback_flips_transport_sign": bool(
                        sign(gf) != 0
                        and sign(gt) != 0
                        and sign(gf) != sign(gt)
                    ),
                    "feedback_dominates_transport_magnitude": bool(
                        abs(gnlf) > abs(gt)
                    ),
                }
            )

        length_rows = [x for x in rows if x["kind"] == "length"]
        diameter_rows = [
            x for x in rows if x["kind"] == "diameter"
        ]
        rec = {
            "drive_scale": float(scale),
            "peak_mV": float(result["soma_mV"][peak_index]),
            "peak_time_ms": float((peak_index + 1) * DT_MS),
            "max_local_depol_mV": float(
                np.max(result["local_depol_mV"])
            ),
            "full_l2": float(np.linalg.norm(full)),
            "frozen_transport_l2": float(
                np.linalg.norm(frozen)
            ),
            "nmda_feedback_l2": float(
                np.linalg.norm(feedback)
            ),
            "full_vs_transport_cosine": cosine(full, frozen),
            "feedback_to_transport_norm_ratio": float(
                np.linalg.norm(feedback)
                / (np.linalg.norm(frozen) + 1e-30)
            ),
            "sign_flips_full_vs_transport": int(
                sum(
                    x["feedback_flips_transport_sign"]
                    for x in rows
                )
            ),
            "length_sign_flips_full_vs_transport": int(
                sum(
                    x["feedback_flips_transport_sign"]
                    for x in length_rows
                )
            ),
            "diameter_sign_flips_full_vs_transport": int(
                sum(
                    x["feedback_flips_transport_sign"]
                    for x in diameter_rows
                )
            ),
            "diameter_full_positive": int(
                sum(x["full_sign"] > 0 for x in diameter_rows)
            ),
            "diameter_transport_positive": int(
                sum(x["frozen_sign"] > 0 for x in diameter_rows)
            ),
            "length_full_positive": int(
                sum(x["full_sign"] > 0 for x in length_rows)
            ),
            "length_transport_positive": int(
                sum(x["frozen_sign"] > 0 for x in length_rows)
            ),
            "feedback_dominates_count": int(
                sum(
                    x["feedback_dominates_transport_magnitude"]
                    for x in rows
                )
            ),
            "rows": rows,
            "all_converged": bool(result["all_converged"]),
            "max_newton_iterations": int(
                result["max_newton_iterations"]
            ),
        }
        records.append(rec)

    all_converged = all(
        x["all_converged"] for x in records
    )
    if reference_error is None:
        raise RuntimeError("1x reference check missing")

    strict = all_converged and reference_error <= 1e-12
    summary = {
        "object": (
            "exact geometry-gradient decomposition into frozen-current "
            "transport and voltage-dependent NMDA feedback"
        ),
        "fci_commit": FCI_COMMIT,
        "cell": {
            "morphology_identifier": model[
                "morphology_identifier"
            ],
            "species": model["species"],
            "layer": model["layer"],
            "compartments": int(len(graph["nodes"])),
            "branch": section_name(branch),
            "branch_directions": int(len(directions)),
        },
        "protocol": {
            "drive_scales": list(DRIVE_SCALES),
            "site_x": list(SITE_X),
            "decomposition": {
                "frozen_current_transport": (
                    "differentiate compiled passive recurrence and "
                    "geometry-dependent input map X while replaying the "
                    "base nonlinear synaptic current waveform"
                ),
                "nmda_feedback": (
                    "full implicit tangent minus frozen-current "
                    "transport tangent"
                ),
            },
        },
        "full_tangent_reference_relative_error_at_1x": float(
            reference_error
        ),
        "records": records,
        "classification": (
            "TRANSPORT_NMDA_GRADIENT_DECOMPOSITION_VALID"
            if strict
            else "TRANSPORT_NMDA_GRADIENT_DECOMPOSITION_FAILED"
        ),
        "stopping_line": (
            "The split is a causal sensitivity decomposition around the "
            "same base trajectory. Frozen-current transport is not a "
            "separate biological neuron; it is the counterfactual where "
            "the nonlinear current waveform is clamped while geometry "
            "changes."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Operaattori transport / NMDA gradient decomposition")
    print(
        f"cell {model['morphology_identifier']} "
        f"branch={section_name(branch)} "
        f"directions={len(directions)}"
    )
    print(
        "full tangent reference error at 1x: "
        f"{reference_error:.3e}"
    )
    for rec in records:
        print(
            f"drive {rec['drive_scale']:>4.2f}: "
            f"local={rec['max_local_depol_mV']:.4g} mV "
            f"||full||={rec['full_l2']:.5g} "
            f"||transport||={rec['frozen_transport_l2']:.5g} "
            f"||feedback||={rec['nmda_feedback_l2']:.5g} "
            f"cos={rec['full_vs_transport_cosine']:+.3f} "
            f"signflips={rec['sign_flips_full_vs_transport']}/"
            f"{len(directions)} "
            f"D+ full/transport="
            f"{rec['diameter_full_positive']}/"
            f"{rec['diameter_transport_positive']}"
        )
    print(f"classification: {summary['classification']}")
    if not strict:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
