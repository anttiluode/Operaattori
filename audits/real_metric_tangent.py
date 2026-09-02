from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "audits"
if str(AUDITS) not in sys.path:
    sys.path.insert(0, str(AUDITS))

from causal_nonlinear_graph import (
    DT_MS,
    REST_MV,
    current_derivative,
    measure_universal_probe_template,
    solve_site_step,
    synapse_law,
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

EPS_LOG = 2e-5
WINDOW_MS = 40.0


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def relerr(actual: np.ndarray, pred: np.ndarray) -> float:
    return rms(np.asarray(pred) - np.asarray(actual)) / (rms(actual) + 1e-30)


def _metric_coefficients(kind: str) -> tuple[float, float]:
    if kind == "length":
        return 1.0, -1.0
    if kind == "diameter":
        return 1.0, 2.0
    if kind == "pose":
        return 0.0, 0.0
    raise ValueError(f"unknown metric kind {kind!r}")


def assemble_metric_graph(
    graph: dict,
    *,
    target_node: int | None = None,
    kind: str = "pose",
    log_scale: float = 0.0,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Reassemble the direct cable graph after one local log-metric change."""
    n = len(graph["nodes"])
    if target_node is not None and not (0 <= int(target_node) < n):
        raise IndexError("target_node outside graph")

    mem_power, axial_power = _metric_coefficients(kind)
    mem_scale = np.ones(n, dtype=float)
    axial_scale = np.ones(n, dtype=float)
    if target_node is not None:
        q = math.exp(float(log_scale))
        mem_scale[int(target_node)] = q ** mem_power
        axial_scale[int(target_node)] = q ** axial_power

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []

    for node in graph["nodes"]:
        i = int(node["node"])
        rows.append(i)
        cols.append(i)
        vals.append(float(node["leak_uS"]) * mem_scale[i])

    for junction in graph["junctions"]:
        edges = [
            (
                int(edge["node"]),
                float(edge["half_axial_uS"])
                * axial_scale[int(edge["node"])],
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
        (vals, (rows, cols)), shape=(n, n), dtype=float
    ).tocsr()
    G.sum_duplicates()

    C = np.asarray(
        [
            float(node["c_nF"]) * mem_scale[int(node["node"])]
            for node in graph["nodes"]
        ],
        dtype=float,
    )
    return G, C


def metric_tangent(
    graph: dict,
    target_node: int,
    kind: str,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Exact dG/dlog(q), dC/dlog(q) for one local metric parameter."""
    n = len(graph["nodes"])
    target = int(target_node)
    if not (0 <= target < n):
        raise IndexError("target_node outside graph")
    mem_power, axial_power = _metric_coefficients(kind)

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []

    if mem_power:
        node = graph["nodes"][target]
        rows.append(target)
        cols.append(target)
        vals.append(mem_power * float(node["leak_uS"]))

    for junction in graph["junctions"]:
        nodes = np.asarray(
            [int(edge["node"]) for edge in junction["edges"]],
            dtype=int,
        )
        g = np.asarray(
            [float(edge["half_axial_uS"]) for edge in junction["edges"]],
            dtype=float,
        )
        dg = np.where(nodes == target, axial_power * g, 0.0)
        anchor = junction["anchored_node"]

        if anchor is not None:
            anchor = int(anchor)
            for node, dgi in zip(nodes, dg):
                node = int(node)
                if node == anchor or dgi == 0.0:
                    continue
                rows.extend([node, anchor, node, anchor])
                cols.extend([node, anchor, anchor, node])
                vals.extend([dgi, dgi, -dgi, -dgi])
            continue

        if len(g) <= 1:
            continue

        total_g = float(np.sum(g))
        dtotal = float(np.sum(dg))
        for i, (ni, gi, dgi) in enumerate(zip(nodes, g, dg)):
            diag = (
                dgi
                - 2.0 * gi * dgi / total_g
                + gi * gi * dtotal / (total_g * total_g)
            )
            if diag != 0.0:
                rows.append(int(ni))
                cols.append(int(ni))
                vals.append(float(diag))

            for j, (nj, gj, dgj) in enumerate(zip(nodes, g, dg)):
                if i == j:
                    continue
                off = (
                    -(dgi * gj + gi * dgj) / total_g
                    + gi * gj * dtotal / (total_g * total_g)
                )
                if off != 0.0:
                    rows.append(int(ni))
                    cols.append(int(nj))
                    vals.append(float(off))

    dG = sparse.coo_matrix(
        (vals, (rows, cols)), shape=(n, n), dtype=float
    ).tocsr()
    dG.sum_duplicates()

    dC = np.zeros(n, dtype=float)
    if mem_power:
        dC[target] = mem_power * float(graph["nodes"][target]["c_nF"])
    return dG, dC


def compile_sparse(
    G: sparse.spmatrix,
    C: np.ndarray,
    site_nodes: np.ndarray,
    soma_node: int,
) -> dict:
    C = np.asarray(C, dtype=float)
    n = len(C)
    d = C / DT_MS
    A = (
        G
        + sparse.diags(d, offsets=0, shape=(n, n), format="csr")
    ).tocsc()
    lu = splu(A)

    B = np.zeros((n, len(site_nodes)), dtype=float)
    for j, node in enumerate(site_nodes):
        B[int(node), int(j)] = 1.0
    X = lu.solve(B)
    R = X[site_nodes, :]
    return {
        "d": d,
        "A": A,
        "lu": lu,
        "B": B,
        "X": X,
        "R": R,
        "site_nodes": np.asarray(site_nodes, dtype=int),
        "soma_node": int(soma_node),
    }


def simulate(
    G: sparse.spmatrix,
    C: np.ndarray,
    site_nodes: np.ndarray,
    soma_node: int,
    ga: np.ndarray,
    gn: np.ndarray,
) -> dict:
    comp = compile_sparse(G, C, site_nodes, soma_node)
    d = comp["d"]
    lu = comp["lu"]
    X = comp["X"]
    R = comp["R"]

    state = np.zeros(len(C), dtype=float)
    previous_z = np.zeros(len(site_nodes), dtype=float)
    soma = np.zeros(ga.shape[1], dtype=float)
    local = np.zeros((len(site_nodes), ga.shape[1]), dtype=float)
    max_iter = 0
    max_consistency = 0.0
    all_converged = True

    for ti in range(ga.shape[1]):
        passive = lu.solve(d * state)
        solved = solve_site_step(
            passive[site_nodes],
            R,
            ga[:, ti],
            gn[:, ti],
            previous_z,
        )
        max_iter = max(max_iter, int(solved["iterations"]))
        all_converged = all_converged and bool(solved["converged"])
        z = np.asarray(solved["z"], dtype=float)
        current = synapse_law(
            REST_MV + z,
            ga[:, ti],
            gn[:, ti],
        )
        state = passive + X @ current
        max_consistency = max(
            max_consistency,
            float(np.max(np.abs(state[site_nodes] - z))),
        )
        soma[ti] = state[soma_node]
        local[:, ti] = state[site_nodes]
        previous_z = state[site_nodes].copy()

    return {
        "soma_mV": soma,
        "local_depol_mV": local,
        "all_converged": bool(all_converged),
        "max_newton_iterations": int(max_iter),
        "max_site_consistency_mV": float(max_consistency),
    }


def simulate_with_metric_tangents(
    G: sparse.spmatrix,
    C: np.ndarray,
    site_nodes: np.ndarray,
    soma_node: int,
    ga: np.ndarray,
    gn: np.ndarray,
    directions: list[tuple[sparse.csr_matrix, np.ndarray]],
) -> dict:
    comp = compile_sparse(G, C, site_nodes, soma_node)
    d = comp["d"]
    lu = comp["lu"]
    X = comp["X"]
    R = comp["R"]

    n = len(C)
    k = len(directions)
    dA = []
    dd = []
    dX = []
    dR = []
    for dG, dC in directions:
        ddi = np.asarray(dC, dtype=float) / DT_MS
        dAi = (
            dG
            + sparse.diags(ddi, offsets=0, shape=(n, n), format="csr")
        ).tocsr()
        dXi = lu.solve(-(dAi @ X))
        dA.append(dAi)
        dd.append(ddi)
        dX.append(dXi)
        dR.append(dXi[site_nodes, :])

    state = np.zeros(n, dtype=float)
    dstate = np.zeros((n, k), dtype=float)
    previous_z = np.zeros(len(site_nodes), dtype=float)

    soma = np.zeros(ga.shape[1], dtype=float)
    soma_tangent = np.zeros((k, ga.shape[1]), dtype=float)
    max_iter = 0
    max_consistency = 0.0
    max_tangent_consistency = 0.0
    all_converged = True

    for ti in range(ga.shape[1]):
        old_state = state
        old_dstate = dstate

        passive = lu.solve(d * old_state)
        tangent_rhs = np.empty((n, k), dtype=float)
        for pi in range(k):
            tangent_rhs[:, pi] = (
                dd[pi] * old_state
                + d * old_dstate[:, pi]
                - dA[pi] @ passive
            )
        dpassive = lu.solve(tangent_rhs)

        solved = solve_site_step(
            passive[site_nodes],
            R,
            ga[:, ti],
            gn[:, ti],
            previous_z,
        )
        max_iter = max(max_iter, int(solved["iterations"]))
        all_converged = all_converged and bool(solved["converged"])

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
                dpassive[site_nodes, pi]
                + dR[pi] @ current
            )
        dz = np.linalg.solve(K, dz_rhs)
        dcurrent = jv[:, None] * dz

        state = passive + X @ current
        dstate = np.empty_like(old_dstate)
        for pi in range(k):
            dstate[:, pi] = (
                dpassive[:, pi]
                + dX[pi] @ current
                + X @ dcurrent[:, pi]
            )

        max_consistency = max(
            max_consistency,
            float(np.max(np.abs(state[site_nodes] - z))),
        )
        max_tangent_consistency = max(
            max_tangent_consistency,
            float(
                np.max(
                    np.abs(dstate[site_nodes, :].T - dz.T)
                )
            ),
        )

        soma[ti] = state[soma_node]
        soma_tangent[:, ti] = dstate[soma_node, :]
        previous_z = state[site_nodes].copy()

    return {
        "soma_mV": soma,
        "soma_tangent_mV_per_logscale": soma_tangent,
        "all_converged": bool(all_converged),
        "max_newton_iterations": int(max_iter),
        "max_site_consistency_mV": float(max_consistency),
        "max_site_tangent_consistency_mV_per_logscale": float(
            max_tangent_consistency
        ),
    }


def sparse_relerr(actual: sparse.spmatrix, pred: sparse.spmatrix) -> float:
    diff = (pred - actual).tocsr()
    return float(
        np.sqrt(np.sum(diff.data * diff.data))
        / (
            np.sqrt(np.sum(actual.data * actual.data))
            + 1e-30
        )
    )


def find_cell1125(panel: list[dict]) -> dict:
    hits = []
    for model in panel:
        haystack = " ".join(str(v) for v in model.values())
        if "1125" in haystack:
            hits.append(model)
    if len(hits) != 1:
        raise RuntimeError(f"expected one cell-1125 panel row, found {len(hits)}")
    return hits[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/real_metric_tangent.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source not pinned")

    setup_neuron(fci_root)
    panel = load_panel(fci_root)
    model = find_cell1125(panel)
    cell = instantiate_matched_passive(fci_root, model)
    graph = build_compartment_graph(cell)

    rebuilt_G, rebuilt_C = assemble_metric_graph(graph)
    base_graph_G_error = sparse_relerr(graph["G_uS"], rebuilt_G)
    base_graph_C_error = relerr(graph["C_nF"], rebuilt_C)

    branches = choose_branches(cell)[:3]
    branch = branches[0]
    site_nodes = np.asarray(
        [
            node_for_section_x(graph, branch, float(x))
            for x in SITE_X
        ],
        dtype=int,
    )
    soma_node = node_for_section_x(graph, cell.soma[0], 0.5)

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
    ga = ga[:, :ntime]
    gn = gn[:, :ntime]

    directions = []
    labels = []
    graph_checks = []
    for source_index, node in enumerate(site_nodes):
        for kind in ("length", "diameter"):
            dG, dC = metric_tangent(graph, int(node), kind)
            directions.append((dG, dC))
            labels.append(
                {
                    "source_index": int(source_index),
                    "site_x": float(SITE_X[source_index]),
                    "node": int(node),
                    "kind": kind,
                }
            )

            Gp, Cp = assemble_metric_graph(
                graph,
                target_node=int(node),
                kind=kind,
                log_scale=EPS_LOG,
            )
            Gm, Cm = assemble_metric_graph(
                graph,
                target_node=int(node),
                kind=kind,
                log_scale=-EPS_LOG,
            )
            fdG = (Gp - Gm) * (1.0 / (2.0 * EPS_LOG))
            fdC = (Cp - Cm) / (2.0 * EPS_LOG)
            graph_checks.append(
                {
                    **labels[-1],
                    "dG_relative_error": sparse_relerr(fdG, dG),
                    "dC_relative_error": relerr(fdC, dC),
                }
            )

    pose_dG, pose_dC = metric_tangent(
        graph, int(site_nodes[1]), "pose"
    )
    pose_zero = (
        pose_dG.nnz == 0
        and float(np.max(np.abs(pose_dC))) == 0.0
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

    base = simulate(
        graph["G_uS"],
        graph["C_nF"],
        site_nodes,
        soma_node,
        ga,
        gn,
    )
    base_trace_match = relerr(
        base["soma_mV"],
        tangent["soma_mV"],
    )

    trace_checks = []
    for pi, label in enumerate(labels):
        Gp, Cp = assemble_metric_graph(
            graph,
            target_node=int(label["node"]),
            kind=label["kind"],
            log_scale=EPS_LOG,
        )
        Gm, Cm = assemble_metric_graph(
            graph,
            target_node=int(label["node"]),
            kind=label["kind"],
            log_scale=-EPS_LOG,
        )
        plus = simulate(
            Gp, Cp, site_nodes, soma_node, ga, gn
        )
        minus = simulate(
            Gm, Cm, site_nodes, soma_node, ga, gn
        )
        numeric = (
            plus["soma_mV"] - minus["soma_mV"]
        ) / (2.0 * EPS_LOG)
        analytic = tangent[
            "soma_tangent_mV_per_logscale"
        ][pi]
        peak_index = int(np.argmax(base["soma_mV"]))
        trace_checks.append(
            {
                **label,
                "trace_tangent_relative_error": relerr(
                    numeric, analytic
                ),
                "peak_time_ms": float(
                    (peak_index + 1) * DT_MS
                ),
                "base_peak_mV": float(
                    base["soma_mV"][peak_index]
                ),
                "analytic_dpeak_mV_per_logscale": float(
                    analytic[peak_index]
                ),
                "numeric_dpeak_mV_per_logscale": float(
                    numeric[peak_index]
                ),
                "one_percent_peak_change_mV_linearized": float(
                    0.01 * analytic[peak_index]
                ),
                "plus_converged": bool(
                    plus["all_converged"]
                ),
                "minus_converged": bool(
                    minus["all_converged"]
                ),
            }
        )

    # The six finite-difference directions above validate the tangent.
    # Now use that analytic object to map every compartment on the tested
    # branch without doing another finite-difference sweep.
    branch_name = section_name(branch)
    branch_nodes = [
        int(x) for x in graph["node_by_section"][branch_name]
    ]
    branch_directions = []
    branch_meta = []
    for node_index in branch_nodes:
        node = graph["nodes"][node_index]
        for kind in ("length", "diameter"):
            branch_directions.append(
                metric_tangent(graph, node_index, kind)
            )
            branch_meta.append(
                {
                    "node": int(node_index),
                    "seg_index": int(node["seg_index"]),
                    "x": float(node["x"]),
                    "length_um": float(node["length_um"]),
                    "diam_um": float(node["diam_um"]),
                    "kind": kind,
                }
            )

    branch_tangent = simulate_with_metric_tangents(
        graph["G_uS"],
        graph["C_nF"],
        site_nodes,
        soma_node,
        ga,
        gn,
        branch_directions,
    )
    peak_index = int(np.argmax(base["soma_mV"]))
    base_peak = float(base["soma_mV"][peak_index])
    branch_map = []
    for pi, meta in enumerate(branch_meta):
        dpeak = float(
            branch_tangent[
                "soma_tangent_mV_per_logscale"
            ][pi, peak_index]
        )
        branch_map.append(
            {
                **meta,
                "dpeak_mV_per_logscale": dpeak,
                "one_percent_peak_change_mV_linearized": float(
                    0.01 * dpeak
                ),
                "one_percent_peak_change_percent_of_peak": float(
                    100.0 * (0.01 * dpeak) / (base_peak + 1e-30)
                ),
            }
        )

    def strongest(rows, kind):
        subset = [x for x in rows if x["kind"] == kind]
        return max(
            subset,
            key=lambda x: abs(x["dpeak_mV_per_logscale"]),
        )

    length_rows = [x for x in branch_map if x["kind"] == "length"]
    diameter_rows = [x for x in branch_map if x["kind"] == "diameter"]
    branch_summary = {
        "compartments": int(len(branch_nodes)),
        "directions": int(len(branch_map)),
        "peak_time_ms": float((peak_index + 1) * DT_MS),
        "base_peak_mV": base_peak,
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
        "strongest_length": strongest(branch_map, "length"),
        "strongest_diameter": strongest(branch_map, "diameter"),
        "all_converged": bool(branch_tangent["all_converged"]),
        "max_site_tangent_consistency_mV_per_logscale": float(
            branch_tangent[
                "max_site_tangent_consistency_mV_per_logscale"
            ]
        ),
    }

    max_graph_G_error = max(
        x["dG_relative_error"] for x in graph_checks
    )
    max_graph_C_error = max(
        x["dC_relative_error"] for x in graph_checks
    )
    max_trace_error = max(
        x["trace_tangent_relative_error"]
        for x in trace_checks
    )

    strict = (
        base_graph_G_error <= 1e-12
        and base_graph_C_error <= 1e-12
        and max_graph_G_error <= 1e-6
        and max_graph_C_error <= 1e-6
        and max_trace_error <= 2e-4
        and base_trace_match <= 1e-12
        and pose_zero
        and tangent["all_converged"]
        and tangent[
            "max_site_tangent_consistency_mV_per_logscale"
        ] <= 1e-7
        and branch_tangent["all_converged"]
        and branch_tangent[
            "max_site_tangent_consistency_mV_per_logscale"
        ] <= 1e-7
    )
    classification = (
        "REAL_MORPHOLOGY_METRIC_TANGENT_VALID"
        if strict
        else "REAL_MORPHOLOGY_METRIC_TANGENT_FAILED"
    )

    summary = {
        "object": (
            "exact local length/diameter derivatives through the "
            "cell-1125 direct cable compiler and causal NMDA soma trace"
        ),
        "fci_commit": FCI_COMMIT,
        "cell": {
            "order": int(model["order"]),
            "species": model["species"],
            "layer": model["layer"],
            "morphology_identifier": model["morphology_identifier"],
            "compartments": int(len(graph["nodes"])),
            "branch": section_name(branch),
            "site_nodes": [int(x) for x in site_nodes],
            "site_x": list(SITE_X),
        },
        "protocol": {
            "parameterization": "theta = log(local metric scale)",
            "directions": (
                "local segment length and diameter at each of the "
                "three branch sites; pose is exact-zero control"
            ),
            "finite_difference_log_step": EPS_LOG,
            "dt_ms": DT_MS,
            "window_ms": WINDOW_MS,
            "timing_ms": [0.0, 5.0, 10.0],
            "conductance_program": "released HUMAN_PROBE template",
            "topology_fixed": True,
            "no_neuron_rerun_for_metric_perturbations": True,
        },
        "base_reassembly": {
            "G_relative_error": base_graph_G_error,
            "C_relative_error": base_graph_C_error,
        },
        "graph_tangent_checks": graph_checks,
        "causal_tangent": {
            "base_trace_match_relative_error": base_trace_match,
            "all_converged": tangent["all_converged"],
            "max_newton_iterations": tangent[
                "max_newton_iterations"
            ],
            "max_site_consistency_mV": tangent[
                "max_site_consistency_mV"
            ],
            "max_site_tangent_consistency_mV_per_logscale": tangent[
                "max_site_tangent_consistency_mV_per_logscale"
            ],
            "trace_checks": trace_checks,
        },
        "pose_zero_control": bool(pose_zero),
        "branch_sensitivity_map": {
            "summary": branch_summary,
            "rows": branch_map,
        },
        "aggregate": {
            "max_dG_relative_error": float(
                max_graph_G_error
            ),
            "max_dC_relative_error": float(
                max_graph_C_error
            ),
            "max_soma_trace_tangent_relative_error": float(
                max_trace_error
            ),
        },
        "classification": classification,
        "stopping_line": (
            "This validates differentiation of the compiled cable abstraction. "
            "It does not establish that biological dendrites optimize soma peak, "
            "nor does it model intracellular organelles unless they change the "
            "effective membrane, axial, channel, or topology state."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Operaattori real morphology metric tangent")
    print(
        f"cell: {model['morphology_identifier']} "
        f"N={len(graph['nodes'])}"
    )
    print(
        "base G/C reassembly error: "
        f"{base_graph_G_error:.3e} / {base_graph_C_error:.3e}"
    )
    print(
        "max graph dG/dC tangent error: "
        f"{max_graph_G_error:.3e} / {max_graph_C_error:.3e}"
    )
    print(
        "max nonlinear soma tangent error: "
        f"{max_trace_error:.3e}"
    )
    print(
        "pose tangent exact zero: "
        f"{pose_zero}"
    )
    print(
        "branch sensitivity directions: "
        f"{branch_summary['directions']}"
    )
    print(
        "strongest length dpeak: "
        f"{branch_summary['strongest_length']['dpeak_mV_per_logscale']:.6g} "
        f"at x={branch_summary['strongest_length']['x']:.3f}"
    )
    print(
        "strongest diameter dpeak: "
        f"{branch_summary['strongest_diameter']['dpeak_mV_per_logscale']:.6g} "
        f"at x={branch_summary['strongest_diameter']['x']:.3f}"
    )
    print(
        "classification: "
        f"{classification}"
    )

    if not strict:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
