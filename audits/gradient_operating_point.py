from __future__ import annotations

import argparse
import json
import math
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
    measure_universal_probe_template,
    timed_conductances,
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
from real_metric_tangent import (
    find_cell1125,
    metric_tangent,
    simulate,
    simulate_with_metric_tangents,
)

WINDOW_MS = 40.0
DRIVE_SCALES = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
L2_BUDGETS = (0.01, 0.05, 0.10, 0.20)
LOCAL_LOG_STEPS = (0.01, 0.05, 0.10, 0.20)


def assemble_metric_field(
    graph: dict,
    log_length: np.ndarray,
    log_diameter: np.ndarray,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Reassemble G,C for simultaneous fixed-topology metric changes.

    For each compartment i:
        membrane area / leak / capacitance ~ exp(ell_i + rho_i)
        half axial conductance             ~ exp(-ell_i + 2 rho_i)

    where ell_i = d log(length_i) and rho_i = d log(diameter_i).
    """
    n = len(graph["nodes"])
    ell = np.asarray(log_length, dtype=float)
    rho = np.asarray(log_diameter, dtype=float)
    if ell.shape != (n,) or rho.shape != (n,):
        raise ValueError("metric fields must have one value per compartment")

    membrane_scale = np.exp(ell + rho)
    axial_scale = np.exp(-ell + 2.0 * rho)

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []

    for node in graph["nodes"]:
        i = int(node["node"])
        rows.append(i)
        cols.append(i)
        vals.append(
            float(node["leak_uS"]) * float(membrane_scale[i])
        )

    for junction in graph["junctions"]:
        edges = [
            (
                int(edge["node"]),
                float(edge["half_axial_uS"])
                * float(axial_scale[int(edge["node"])]),
            )
            for edge in junction["edges"]
        ]
        anchor = junction["anchored_node"]

        if anchor is not None:
            anchor = int(anchor)
            for node, g in edges:
                if node == anchor:
                    continue
                rows.extend([node, anchor, node, anchor])
                cols.extend([node, anchor, anchor, node])
                vals.extend([g, g, -g, -g])
            continue

        if len(edges) <= 1:
            continue

        total_g = float(sum(g for _, g in edges))
        for i, (ni, gi) in enumerate(edges):
            rows.append(ni)
            cols.append(ni)
            vals.append(gi - gi * gi / total_g)
            for j, (nj, gj) in enumerate(edges):
                if i == j:
                    continue
                rows.append(ni)
                cols.append(nj)
                vals.append(-gi * gj / total_g)

    G = sparse.coo_matrix(
        (vals, (rows, cols)),
        shape=(n, n),
        dtype=float,
    ).tocsr()
    G.sum_duplicates()

    C = np.asarray(
        [
            float(node["c_nF"])
            * float(membrane_scale[int(node["node"])])
            for node in graph["nodes"]
        ],
        dtype=float,
    )
    return G, C


def branch_directions(
    graph: dict,
    branch,
) -> tuple[list[tuple[sparse.csr_matrix, np.ndarray]], list[dict]]:
    name = section_name(branch)
    nodes = [int(x) for x in graph["node_by_section"][name]]
    directions = []
    meta = []
    for node_index in nodes:
        node = graph["nodes"][node_index]
        for kind in ("length", "diameter"):
            directions.append(metric_tangent(graph, node_index, kind))
            meta.append(
                {
                    "node": int(node_index),
                    "seg_index": int(node["seg_index"]),
                    "x": float(node["x"]),
                    "kind": kind,
                    "length_um": float(node["length_um"]),
                    "diam_um": float(node["diam_um"]),
                }
            )
    return directions, meta


def conductance_program(scale: float) -> tuple[np.ndarray, np.ndarray]:
    template = measure_universal_probe_template()
    ga, gn = timed_conductances(
        template["g_ampa_uS"],
        template["g_nmda_raw_uS"],
        (0.0, 5.0, 10.0),
    )
    ntime = min(
        ga.shape[1],
        int(round(WINDOW_MS / DT_MS)),
    )
    return (
        float(scale) * ga[:, :ntime],
        float(scale) * gn[:, :ntime],
    )


def operating_point(
    graph: dict,
    site_nodes: np.ndarray,
    soma_node: int,
    directions,
    meta,
    scale: float,
) -> dict:
    ga, gn = conductance_program(scale)
    base = simulate(
        graph["G_uS"],
        graph["C_nF"],
        site_nodes,
        soma_node,
        ga,
        gn,
    )
    tangent = simulate_with_metric_tangents(
        graph["G_uS"],
        graph["C_nF"],
        site_nodes,
        soma_node,
        ga,
        gn,
        directions,
    )
    peak_index = int(np.argmax(base["soma_mV"]))
    peak_mV = float(base["soma_mV"][peak_index])
    gradient = np.asarray(
        tangent["soma_tangent_mV_per_logscale"][:, peak_index],
        dtype=float,
    )

    rows = []
    for mi, gi in zip(meta, gradient):
        rows.append({**mi, "dpeak_mV_per_logscale": float(gi)})

    length_rows = [x for x in rows if x["kind"] == "length"]
    diameter_rows = [x for x in rows if x["kind"] == "diameter"]
    positive_length = [x for x in length_rows if x["dpeak_mV_per_logscale"] > 0]
    negative_length = [x for x in length_rows if x["dpeak_mV_per_logscale"] < 0]

    strongest = lambda subset: max(
        subset,
        key=lambda x: abs(x["dpeak_mV_per_logscale"]),
    )
    first_negative = (
        min(negative_length, key=lambda x: x["x"])
        if negative_length
        else None
    )

    return {
        "drive_scale": float(scale),
        "peak_mV": peak_mV,
        "peak_time_ms": float((peak_index + 1) * DT_MS),
        "max_local_depol_mV": float(
            np.max(base["local_depol_mV"])
        ),
        "gradient_l2_mV_per_logscale": float(
            np.linalg.norm(gradient)
        ),
        "gradient_linf_mV_per_logscale": float(
            np.max(np.abs(gradient))
        ),
        "length_positive": int(
            sum(x["dpeak_mV_per_logscale"] > 0 for x in length_rows)
        ),
        "length_negative": int(
            sum(x["dpeak_mV_per_logscale"] < 0 for x in length_rows)
        ),
        "diameter_positive": int(
            sum(x["dpeak_mV_per_logscale"] > 0 for x in diameter_rows)
        ),
        "diameter_negative": int(
            sum(x["dpeak_mV_per_logscale"] < 0 for x in diameter_rows)
        ),
        "first_negative_length_x": (
            None if first_negative is None else float(first_negative["x"])
        ),
        "strongest_length": strongest(length_rows),
        "strongest_diameter": strongest(diameter_rows),
        "strongest_positive_length": (
            None
            if not positive_length
            else max(
                positive_length,
                key=lambda x: x["dpeak_mV_per_logscale"],
            )
        ),
        "strongest_negative_length": (
            None
            if not negative_length
            else min(
                negative_length,
                key=lambda x: x["dpeak_mV_per_logscale"],
            )
        ),
        "all_converged": bool(
            base["all_converged"] and tangent["all_converged"]
        ),
        "max_newton_iterations": int(
            max(
                base["max_newton_iterations"],
                tangent["max_newton_iterations"],
            )
        ),
        "gradient": [float(x) for x in gradient],
        "rows": rows,
    }


def vector_to_fields(
    graph: dict,
    meta: list[dict],
    theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(graph["nodes"])
    ell = np.zeros(n, dtype=float)
    rho = np.zeros(n, dtype=float)
    for m, value in zip(meta, np.asarray(theta, dtype=float)):
        if m["kind"] == "length":
            ell[int(m["node"])] += float(value)
        else:
            rho[int(m["node"])] += float(value)
    return ell, rho


def finite_budget_sweep(
    graph: dict,
    site_nodes: np.ndarray,
    soma_node: int,
    meta: list[dict],
    baseline_point: dict,
) -> dict:
    scale = float(baseline_point["drive_scale"])
    ga, gn = conductance_program(scale)
    base = simulate(
        graph["G_uS"],
        graph["C_nF"],
        site_nodes,
        soma_node,
        ga,
        gn,
    )
    peak_index = int(np.argmax(base["soma_mV"]))
    base_fixed = float(base["soma_mV"][peak_index])
    base_peak = float(np.max(base["soma_mV"]))
    g = np.asarray(baseline_point["gradient"], dtype=float)
    norm = float(np.linalg.norm(g))
    unit = g / (norm + 1e-30)

    l2_rows = []
    for budget in L2_BUDGETS:
        for direction_name, sign in (
            ("ascent", 1.0),
            ("descent", -1.0),
        ):
            theta = sign * float(budget) * unit
            ell, rho = vector_to_fields(graph, meta, theta)
            Gq, Cq = assemble_metric_field(graph, ell, rho)
            out = simulate(Gq, Cq, site_nodes, soma_node, ga, gn)
            actual_fixed = (
                float(out["soma_mV"][peak_index]) - base_fixed
            )
            actual_peak = float(np.max(out["soma_mV"])) - base_peak
            predicted = float(np.dot(g, theta))
            max_abs_theta = float(np.max(np.abs(theta)))
            l2_rows.append(
                {
                    "direction": direction_name,
                    "l2_log_metric_budget": float(budget),
                    "predicted_fixed_time_change_mV": predicted,
                    "actual_fixed_time_change_mV": actual_fixed,
                    "actual_peak_change_mV": actual_peak,
                    "actual_peak_change_percent": float(
                        100.0 * actual_peak / (base_peak + 1e-30)
                    ),
                    "prediction_ratio_fixed_time": float(
                        actual_fixed / (predicted + 1e-30)
                    ),
                    "max_local_abs_log_change": max_abs_theta,
                    "max_local_abs_percent_change_approx": float(
                        100.0 * (math.exp(max_abs_theta) - 1.0)
                    ),
                    "all_converged": bool(out["all_converged"]),
                }
            )

    length_indices = [
        i for i, m in enumerate(meta) if m["kind"] == "length"
    ]
    positive_length = [
        i for i in length_indices if g[i] > 0.0
    ]
    negative_length = [
        i for i in length_indices if g[i] < 0.0
    ]
    warm_index = max(
        positive_length,
        key=lambda i: g[i],
    )
    cool_index = min(
        negative_length,
        key=lambda i: g[i],
    )

    local_rows = []
    for label, index in (
        ("lengthen_pre_site_warm", warm_index),
        ("lengthen_post_site_cool", cool_index),
    ):
        for step in LOCAL_LOG_STEPS:
            theta = np.zeros_like(g)
            theta[index] = float(step)
            ell, rho = vector_to_fields(graph, meta, theta)
            Gq, Cq = assemble_metric_field(graph, ell, rho)
            out = simulate(Gq, Cq, site_nodes, soma_node, ga, gn)
            predicted = float(g[index] * step)
            actual_fixed = (
                float(out["soma_mV"][peak_index]) - base_fixed
            )
            actual_peak = float(np.max(out["soma_mV"])) - base_peak
            local_rows.append(
                {
                    "experiment": label,
                    "node": int(meta[index]["node"]),
                    "x": float(meta[index]["x"]),
                    "log_length_step": float(step),
                    "length_percent_change": float(
                        100.0 * (math.exp(step) - 1.0)
                    ),
                    "gradient_mV_per_logscale": float(g[index]),
                    "predicted_fixed_time_change_mV": predicted,
                    "actual_fixed_time_change_mV": actual_fixed,
                    "actual_peak_change_mV": actual_peak,
                    "actual_peak_change_percent": float(
                        100.0 * actual_peak / (base_peak + 1e-30)
                    ),
                    "prediction_ratio_fixed_time": float(
                        actual_fixed / (predicted + 1e-30)
                    ),
                    "all_converged": bool(out["all_converged"]),
                }
            )

    return {
        "drive_scale": scale,
        "base_peak_mV": base_peak,
        "base_peak_time_ms": float((peak_index + 1) * DT_MS),
        "gradient_l2_mV_per_logscale": norm,
        "l2_gradient_following": l2_rows,
        "local_lengthening": local_rows,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/"
            "gradient_operating_point.json"
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

    points = [
        operating_point(
            graph,
            site_nodes,
            soma_node,
            directions,
            meta,
            scale,
        )
        for scale in DRIVE_SCALES
    ]
    baseline = next(x for x in points if x["drive_scale"] == 1.0)
    budgets = finite_budget_sweep(
        graph,
        site_nodes,
        soma_node,
        meta,
        baseline,
    )

    all_converged = (
        all(x["all_converged"] for x in points)
        and all(
            x["all_converged"]
            for x in budgets["l2_gradient_following"]
        )
        and all(
            x["all_converged"]
            for x in budgets["local_lengthening"]
        )
    )

    summary = {
        "object": (
            "operating-point stability and finite-budget usefulness "
            "of the verified cell-1125 branch metric gradient"
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
            "branch_compartments": int(
                len(graph["node_by_section"][section_name(branch)])
            ),
            "site_x": list(SITE_X),
        },
        "protocol": {
            "drive_scales": list(DRIVE_SCALES),
            "timing_ms": [0.0, 5.0, 10.0],
            "window_ms": WINDOW_MS,
            "dt_ms": DT_MS,
            "gradient_target": (
                "soma depolarization at the unperturbed peak time"
            ),
            "finite_geometry": (
                "full graph recompile after simultaneous log-length/"
                "log-diameter changes; no tangent used for the actual"
                " finite response"
            ),
            "l2_log_metric_budgets": list(L2_BUDGETS),
            "local_log_length_steps": list(LOCAL_LOG_STEPS),
        },
        "operating_points": points,
        "finite_budget": budgets,
        "all_converged": bool(all_converged),
        "classification": (
            "GRADIENT_OPERATING_POINT_AUDIT_VALID"
            if all_converged
            else "GRADIENT_OPERATING_POINT_AUDIT_FAILED"
        ),
        "stopping_line": (
            "This audit measures sensitivity and finite response in one "
            "fixed-topology branch. It does not establish a biological "
            "growth rule, a useful learning objective, or a behaviorally "
            "meaningful voltage scale."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Operaattori gradient operating-point audit")
    print(
        f"cell {model['morphology_identifier']} "
        f"branch {section_name(branch)} "
        f"Nbranch={summary['cell']['branch_compartments']}"
    )
    for point in points:
        print(
            "drive "
            f"{point['drive_scale']:>4.2f}: "
            f"peak={point['peak_mV']:.6g} mV "
            f"local={point['max_local_depol_mV']:.6g} mV "
            f"L+/-={point['length_positive']}/"
            f"{point['length_negative']} "
            f"D+/-={point['diameter_positive']}/"
            f"{point['diameter_negative']} "
            f"first L- x={point['first_negative_length_x']} "
            f"||g||={point['gradient_l2_mV_per_logscale']:.6g}"
        )

    print("finite L2 gradient following at drive=1:")
    for row in budgets["l2_gradient_following"]:
        if row["direction"] == "ascent":
            print(
                f"  budget={row['l2_log_metric_budget']:.3f} "
                f"maxlocal≈{row['max_local_abs_percent_change_approx']:.3f}% "
                f"dpeak={row['actual_peak_change_mV']:.6g} mV "
                f"({row['actual_peak_change_percent']:.4g}%) "
                f"pred={row['predicted_fixed_time_change_mV']:.6g}"
            )

    print("local lengthening at drive=1:")
    for row in budgets["local_lengthening"]:
        if abs(row["log_length_step"] - 0.20) < 1e-12:
            print(
                f"  {row['experiment']} x={row['x']:.3f}: "
                f"+{row['length_percent_change']:.2f}% length -> "
                f"dpeak={row['actual_peak_change_mV']:.6g} mV "
                f"({row['actual_peak_change_percent']:.4g}%)"
            )

    print(f"classification: {summary['classification']}")
    if not all_converged:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
