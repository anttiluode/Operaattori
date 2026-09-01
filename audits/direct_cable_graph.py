from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "audits"
if str(AUDITS) not in sys.path:
    sys.path.insert(0, str(AUDITS))

from cross_cell_operator import (
    DT_MS,
    FCI_COMMIT,
    IMPULSE_NA,
    N_BRANCHES,
    POST_MS,
    SITE_X,
    V_REST_MV,
    cell_global_geometry,
    choose_branches,
    git_head,
    instantiate_matched_passive,
    load_panel,
    measure_branch_operator,
    run_noinput,
    section_name,
    setup_neuron,
)

RA_OHM_CM = 150.0
RM_OHM_CM2 = 20000.0
CM_UF_CM2 = 1.0


class UnionFind:
    def __init__(self):
        self.parent = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x

    def find(self, x):
        self.add(x)
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def nrmse(actual: np.ndarray, pred: np.ndarray) -> float:
    return rms(np.asarray(pred) - np.asarray(actual)) / (
        rms(actual) + 1e-30
    )


def half_axial_conductance_uS(
    length_um: float,
    diam_um: float,
) -> float:
    half_length_cm = max(float(length_um), 1e-12) * 0.5 * 1e-4
    radius_cm = max(float(diam_um), 1e-9) * 0.5 * 1e-4
    cross_cm2 = math.pi * radius_cm * radius_cm
    resistance_ohm = (
        RA_OHM_CM * half_length_cm / cross_cm2
    )
    return float(1e6 / resistance_ohm)


def parent_connection_x(sec) -> float | None:
    from neuron import h

    try:
        ref = h.SectionRef(sec=sec)
        if not int(ref.has_parent()):
            return None
    except Exception:
        try:
            if sec.parentseg() is None:
                return None
        except Exception:
            return None

    try:
        return float(h.parent_connection(sec=sec))
    except Exception:
        try:
            return float(sec.parentseg().x)
        except Exception:
            return None


def child_orientation(sec) -> float:
    try:
        return float(sec.orientation())
    except Exception:
        return 0.0


def build_compartment_graph(cell) -> dict:
    from neuron import h

    sections = list(cell.all)
    nodes = []
    node_by_section = {}
    boundary_ids = {}
    uf = UnionFind()

    for sec in sections:
        name = section_name(sec)
        segs = list(sec)
        ids = []
        for si, seg in enumerate(segs):
            node_id = len(nodes)
            ids.append(node_id)
            area_um2 = float(h.area(float(seg.x), sec=sec))
            area_cm2 = area_um2 * 1e-8
            c_nF = (
                float(sec.cm) * area_cm2 * 1000.0
            )
            leak_uS = (
                float(sec.g_pas) * area_cm2 * 1e6
            )
            nodes.append(
                {
                    "node": node_id,
                    "section": name,
                    "seg_index": int(si),
                    "x": float(seg.x),
                    "diam_um": float(seg.diam),
                    "length_um": float(sec.L) / max(int(sec.nseg), 1),
                    "c_nF": c_nF,
                    "leak_uS": leak_uS,
                }
            )
        node_by_section[name] = ids

        for k in range(int(sec.nseg) + 1):
            bid = (name, int(k))
            boundary_ids[(name, int(k))] = bid
            uf.add(bid)

    # A child endpoint can connect either to a parent cable boundary or
    # directly to a parent membrane compartment center (notably soma(0.5)).
    anchor_requests = []
    unsupported_connections = []
    boundary_connections = 0
    center_connections = 0

    for sec in sections:
        name = section_name(sec)
        try:
            pseg = sec.parentseg()
        except Exception:
            pseg = None
        if pseg is None:
            continue

        parent = pseg.sec
        pname = section_name(parent)

        orientation = child_orientation(sec)
        child_k = 0 if orientation < 0.5 else int(sec.nseg)
        child_bid = boundary_ids[(name, child_k)]

        px = parent_connection_x(sec)
        if px is None:
            raise RuntimeError(
                f"could not recover parent connection for {name}"
            )

        raw = px * int(parent.nseg)
        parent_k = int(round(raw))
        parent_k = max(0, min(int(parent.nseg), parent_k))
        boundary_x = parent_k / int(parent.nseg)
        boundary_mismatch = abs(px - boundary_x)

        if boundary_mismatch <= 1e-6:
            uf.union(
                child_bid,
                boundary_ids[(pname, parent_k)],
            )
            boundary_connections += 1
            continue

        parent_seg = parent(float(px))
        parent_center_x = float(parent_seg.x)
        center_mismatch = abs(px - parent_center_x)
        if center_mismatch <= 1e-6:
            parent_ids = node_by_section[pname]
            anchor_node = min(
                parent_ids,
                key=lambda idx: abs(
                    nodes[idx]["x"] - parent_center_x
                ),
            )
            anchor_requests.append(
                (child_bid, int(anchor_node))
            )
            center_connections += 1
            continue

        unsupported_connections.append(
            {
                "child": name,
                "parent": pname,
                "parent_x": px,
                "nearest_boundary_x": boundary_x,
                "nearest_segment_center_x": parent_center_x,
                "boundary_mismatch": boundary_mismatch,
                "center_mismatch": center_mismatch,
            }
        )

    if unsupported_connections:
        raise RuntimeError(
            "connection is neither a compartment boundary nor a membrane "
            "center in the locked discretization: "
            + json.dumps(unsupported_connections[:10])
        )

    # Every segment center connects to its two boundaries through half of the
    # segment's axial resistance.
    incident = defaultdict(list)
    for node in nodes:
        name = node["section"]
        si = node["seg_index"]
        g = half_axial_conductance_uS(
            node["length_um"],
            node["diam_um"],
        )
        left = uf.find(boundary_ids[(name, si)])
        right = uf.find(boundary_ids[(name, si + 1)])
        incident[left].append((node["node"], g))
        incident[right].append((node["node"], g))

    # Resolve direct-to-membrane anchors after all boundary unions are known.
    anchored = {}
    for bid, anchor_node in anchor_requests:
        root = uf.find(bid)
        if root in anchored and anchored[root] != anchor_node:
            raise RuntimeError(
                "one zero-junction group was anchored to two membrane centers"
            )
        anchored[root] = anchor_node

    n = len(nodes)
    rows = []
    cols = []
    vals = []

    for node in nodes:
        rows.append(node["node"])
        cols.append(node["node"])
        vals.append(node["leak_uS"])

    junction_degree_counts = defaultdict(int)
    anchored_junctions = 0

    for root, edges in incident.items():
        junction_degree_counts[len(edges)] += 1

        if root in anchored:
            # The zero-length connection is clamped to an existing membrane
            # center. Each attached half-cable therefore contributes a direct
            # axial conductance to that center.
            anchor = anchored[root]
            anchored_junctions += 1
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
        [node["c_nF"] for node in nodes],
        dtype=float,
    )
    if np.any(C <= 0) or not np.all(np.isfinite(C)):
        raise RuntimeError("non-positive/non-finite compartment capacitance")

    return {
        "nodes": nodes,
        "node_by_section": node_by_section,
        "G_uS": G,
        "C_nF": C,
        "boundary_connections": int(boundary_connections),
        "center_connections": int(center_connections),
        "anchored_junctions": int(anchored_junctions),
        "junction_degree_counts": {
            str(k): int(v)
            for k, v in sorted(
                junction_degree_counts.items()
            )
        },
    }


def node_for_section_x(
    graph: dict,
    sec,
    x: float,
) -> int:
    name = section_name(sec)
    ids = graph["node_by_section"][name]

    # Match NEURON's actual Segment selected by sec(x), not raw x.
    seg = sec(float(x))
    target_x = float(seg.x)

    return min(
        ids,
        key=lambda idx: abs(
            graph["nodes"][idx]["x"] - target_x
        ),
    )


def solve_all_sources(
    graph: dict,
    cell,
    branches: list[dict],
    reference_time: np.ndarray,
) -> dict:
    G = graph["G_uS"]
    C = graph["C_nF"]
    n = len(C)

    sources = []
    for branch in branches:
        for source_index, x in enumerate(SITE_X):
            sources.append(
                {
                    "branch_index": int(
                        branch["branch_index"]
                    ),
                    "source_index": int(source_index),
                    "node": node_for_section_x(
                        graph, branch["sec"], x
                    ),
                }
            )
    ns = len(sources)
    if ns != N_BRANCHES * len(SITE_X):
        raise RuntimeError("unexpected source count")

    outputs = {}
    for branch in branches:
        for target_index, x in enumerate(SITE_X):
            outputs[
                (int(branch["branch_index"]), int(target_index))
            ] = node_for_section_x(
                graph, branch["sec"], x
            )
    soma_node = node_for_section_x(
        graph, cell.soma[0], 0.5
    )

    d = C / DT_MS
    A = (
        G
        + sparse.diags(
            d, offsets=0, shape=(n, n), format="csr"
        )
    ).tocsc()
    lu = splu(A)

    reference_time = np.asarray(reference_time, dtype=float)
    if reference_time.ndim != 1 or len(reference_time) == 0:
        raise ValueError("reference_time must be a nonempty vector")

    sample_steps = np.rint(
        reference_time / DT_MS
    ).astype(int)
    if np.any(sample_steps < 0):
        raise RuntimeError("negative graph sample step")
    if not np.allclose(
        reference_time,
        sample_steps * DT_MS,
        rtol=0,
        atol=1e-9,
    ):
        raise RuntimeError(
            "reference times are not on the locked 0.05 ms grid"
        )

    nt = len(reference_time)
    max_step = int(np.max(sample_steps))
    step_to_columns = defaultdict(list)
    for out_index, step in enumerate(sample_steps):
        step_to_columns[int(step)].append(int(out_index))

    state = np.zeros((n, ns), dtype=float)

    soma = np.zeros((ns, nt), dtype=float)
    local = {
        key: np.zeros((ns, nt), dtype=float)
        for key in outputs
    }

    def record(step: int) -> None:
        for out_index in step_to_columns.get(step, []):
            soma[:, out_index] = state[soma_node, :]
            for key, node in outputs.items():
                local[key][:, out_index] = state[node, :]

    record(0)

    input_matrix = np.zeros((n, ns), dtype=float)
    for col, src in enumerate(sources):
        input_matrix[src["node"], col] = IMPULSE_NA

    for step in range(1, max_step + 1):
        rhs = d[:, None] * state
        if step == 1:
            rhs += input_matrix
        state = lu.solve(rhs)
        record(step)

    return {
        "sources": sources,
        "soma": soma,
        "local": local,
        "time": reference_time.copy(),
    }


def make_branch_rows(cell) -> list[dict]:
    rows = []
    for bi, sec in enumerate(choose_branches(cell)):
        rows.append(
            {
                "branch_index": int(bi),
                "section": section_name(sec),
                "sec": sec,
            }
        )
    return rows


def graph_operator_for_branch(
    solved: dict,
    branch_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    nt = len(solved["time"])
    G = np.zeros((3, 3, nt), dtype=float)
    T = np.zeros((3, nt), dtype=float)

    source_cols = {}
    for col, src in enumerate(solved["sources"]):
        if int(src["branch_index"]) == int(branch_index):
            source_cols[int(src["source_index"])] = col

    for source_index in range(3):
        col = source_cols[source_index]
        T[source_index] = solved["soma"][col]
        for target_index in range(3):
            G[target_index, source_index] = solved["local"][
                (int(branch_index), int(target_index))
            ][col]

    # NEURON reference divides a 0.001 nA response by 0.001.  Do the same.
    return G / IMPULSE_NA, T / IMPULSE_NA


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/direct_cable_graph.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source not pinned")

    setup_neuron(fci_root)
    panel = load_panel(fci_root)

    branch_scores = []
    cell_receipts = []

    for ci, model in enumerate(panel):
        cell = instantiate_matched_passive(
            fci_root, model
        )
        branches = make_branch_rows(cell)
        noinput = run_noinput(cell, branches)

        graph = build_compartment_graph(cell)
        solved = solve_all_sources(
            graph, cell, branches, noinput["t"]
        )

        cell_scores = []
        for branch in branches:
            actual_G, actual_T = measure_branch_operator(
                cell, branch, noinput
            )
            pred_G, pred_T = graph_operator_for_branch(
                solved, branch["branch_index"]
            )

            eg = nrmse(actual_G, pred_G)
            et = nrmse(actual_T, pred_T)
            joint = 0.5 * (eg + et)
            row = {
                "cell_order": int(model["order"]),
                "species": model["species"],
                "layer": model["layer"],
                "morphology_identifier": model[
                    "morphology_identifier"
                ],
                "branch_index": int(
                    branch["branch_index"]
                ),
                "section": branch["section"],
                "G_nrmse": float(eg),
                "T_nrmse": float(et),
                "joint_nrmse": float(joint),
            }
            branch_scores.append(row)
            cell_scores.append(row)

        cell_joint = float(
            np.median(
                [x["joint_nrmse"] for x in cell_scores]
            )
        )
        cell_receipts.append(
            {
                "cell_order": int(model["order"]),
                "species": model["species"],
                "layer": model["layer"],
                "morphology_identifier": model[
                    "morphology_identifier"
                ],
                "compartments": int(
                    len(graph["nodes"])
                ),
                "boundary_connections": graph[
                    "boundary_connections"
                ],
                "center_connections": graph[
                    "center_connections"
                ],
                "anchored_junctions": graph[
                    "anchored_junctions"
                ],
                "junction_degree_counts": graph[
                    "junction_degree_counts"
                ],
                "median_G_nrmse": float(
                    np.median(
                        [x["G_nrmse"] for x in cell_scores]
                    )
                ),
                "median_T_nrmse": float(
                    np.median(
                        [x["T_nrmse"] for x in cell_scores]
                    )
                ),
                "median_joint_nrmse": cell_joint,
            }
        )

        print(
            f"[{ci+1:02d}/24] "
            f"{model['species']:5s} "
            f"{model['morphology_identifier']:>12s} "
            f"N={len(graph['nodes']):5d} "
            f"joint={cell_joint:.4f}"
        )

    g = np.asarray(
        [x["G_nrmse"] for x in branch_scores],
        dtype=float,
    )
    t = np.asarray(
        [x["T_nrmse"] for x in branch_scores],
        dtype=float,
    )
    j = np.asarray(
        [x["joint_nrmse"] for x in branch_scores],
        dtype=float,
    )
    cell_j = np.asarray(
        [x["median_joint_nrmse"] for x in cell_receipts],
        dtype=float,
    )

    cells_under_10 = int(
        np.sum(cell_j <= 0.10)
    )
    worst = max(
        cell_receipts,
        key=lambda x: x["median_joint_nrmse"],
    )

    median_joint = float(np.median(j))
    median_g = float(np.median(g))
    median_t = float(np.median(t))
    median_cell = float(np.median(cell_j))

    strict = (
        median_joint <= 0.05
        and median_g <= 0.05
        and median_t <= 0.05
        and median_cell <= 0.06
        and cells_under_10 >= 20
    )

    if strict:
        classification = (
            "MORPHOLOGY_GRAPH_GENERATES_PASSIVE_OPERATOR"
        )
    elif median_joint <= 0.10:
        classification = (
            "MORPHOLOGY_GRAPH_APPROXIMATES_PASSIVE_OPERATOR"
        )
    else:
        classification = (
            "HAND_BUILT_CABLE_DISCRETIZATION_INADEQUATE"
        )

    summary = {
        "object": (
            "direct morphology-graph construction of matched-passive local "
            "Green and site-to-soma transport operators"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "cells": 24,
            "branches_per_cell": N_BRANCHES,
            "sites_per_branch": len(SITE_X),
            "site_x": list(SITE_X),
            "Ra_ohm_cm": RA_OHM_CM,
            "Rm_ohm_cm2": RM_OHM_CM2,
            "Cm_uF_cm2": CM_UF_CM2,
            "dt_ms": DT_MS,
            "post_ms": POST_MS,
            "impulse_nA": IMPULSE_NA,
            "no_cross_cell_fitting": True,
            "junction_elimination": (
                "zero-capacitance boundary Schur complement"
            ),
            "time_step": "backward Euler",
            "thresholds_locked_before_run": {
                "median_joint_nrmse_max": 0.05,
                "median_G_nrmse_max": 0.05,
                "median_T_nrmse_max": 0.05,
                "median_cell_joint_nrmse_max": 0.06,
                "cells_under_0p10_min": 20,
                "approximate_classification_joint_max": 0.10,
            },
        },
        "aggregate": {
            "branch_operator_packs": int(
                len(branch_scores)
            ),
            "median_G_nrmse": median_g,
            "median_T_nrmse": median_t,
            "median_joint_nrmse": median_joint,
            "median_cell_joint_nrmse": median_cell,
            "cells_joint_nrmse_le_0p10": cells_under_10,
            "worst_cell": worst,
        },
        "cells": cell_receipts,
        "branches": branch_scores,
        "classification": classification,
        "stopping_line": (
            "Do not align traces or fit gains/cell-specific corrections. "
            "If this fails, diagnose cable discretization before changing "
            "the biology or adding learning."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori direct morphology-graph audit")
    print()
    print(
        "median joint NRMSE:             "
        f"{median_joint:.4f}"
    )
    print(
        "median G NRMSE:                 "
        f"{median_g:.4f}"
    )
    print(
        "median T NRMSE:                 "
        f"{median_t:.4f}"
    )
    print(
        "median cell joint NRMSE:        "
        f"{median_cell:.4f}"
    )
    print(
        "cells <= 0.10:                  "
        f"{cells_under_10} / 24"
    )
    print(
        "worst cell:                     "
        f"{worst['cell_order']} "
        f"{worst['species']} "
        f"{worst['morphology_identifier']} "
        f"{worst['median_joint_nrmse']:.4f}"
    )
    print(f"classification: {classification}")

    assert len(branch_scores) == 144


if __name__ == "__main__":
    main()
