from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "audits"
if str(AUDITS) not in sys.path:
    sys.path.insert(0, str(AUDITS))

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
    operating_point,
)
from real_metric_tangent import (
    find_cell1125,
    metric_tangent,
)

PROBE_X = (0.10, 0.25, 0.50, 0.75)


def probe_directions(graph: dict, branch) -> tuple[list, list[dict]]:
    directions = []
    meta = []
    seen = {}
    for requested_x in PROBE_X:
        node_index = int(
            node_for_section_x(graph, branch, float(requested_x))
        )
        node = graph["nodes"][node_index]
        actual_x = float(node["x"])
        seen.setdefault(node_index, []).append(float(requested_x))
        for kind in ("length", "diameter"):
            directions.append(
                metric_tangent(graph, node_index, kind)
            )
            meta.append(
                {
                    "node": node_index,
                    "seg_index": int(node["seg_index"]),
                    "x": actual_x,
                    "requested_x": float(requested_x),
                    "kind": kind,
                    "length_um": float(node["length_um"]),
                    "diam_um": float(node["diam_um"]),
                }
            )
    duplicate_nodes = {
        str(k): v for k, v in seen.items() if len(v) > 1
    }
    return directions, meta, duplicate_nodes


def sign(v: float, eps: float = 1e-12) -> int:
    if v > eps:
        return 1
    if v < -eps:
        return -1
    return 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/"
            "gradient_branch_panel.json"
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
    branches = choose_branches(cell)

    records = []
    duplicate_probe_nodes = {}
    for bi, branch in enumerate(branches):
        site_nodes = np.asarray(
            [
                node_for_section_x(graph, branch, float(x))
                for x in SITE_X
            ],
            dtype=int,
        )
        directions, meta, duplicates = probe_directions(
            graph, branch
        )
        if duplicates:
            duplicate_probe_nodes[str(bi)] = duplicates

        for scale in DRIVE_SCALES:
            point = operating_point(
                graph,
                site_nodes,
                node_for_section_x(graph, cell.soma[0], 0.5),
                directions,
                meta,
                scale,
            )
            rows = []
            for row in point["rows"]:
                rows.append(
                    {
                        "requested_x": float(row["requested_x"]),
                        "actual_x": float(row["x"]),
                        "node": int(row["node"]),
                        "kind": row["kind"],
                        "dpeak_mV_per_logscale": float(
                            row["dpeak_mV_per_logscale"]
                        ),
                        "sign": sign(
                            row["dpeak_mV_per_logscale"]
                        ),
                    }
                )
            records.append(
                {
                    "branch_index": int(bi),
                    "section": section_name(branch),
                    "drive_scale": float(scale),
                    "peak_mV": float(point["peak_mV"]),
                    "max_local_depol_mV": float(
                        point["max_local_depol_mV"]
                    ),
                    "gradient_l2_mV_per_logscale": float(
                        point["gradient_l2_mV_per_logscale"]
                    ),
                    "rows": rows,
                    "all_converged": bool(point["all_converged"]),
                }
            )

    by_key = {}
    for record in records:
        for row in record["rows"]:
            key = (
                record["branch_index"],
                row["requested_x"],
                row["kind"],
            )
            by_key.setdefault(key, []).append(
                (
                    record["drive_scale"],
                    row["dpeak_mV_per_logscale"],
                )
            )

    trajectories = []
    for key, values in sorted(by_key.items()):
        values.sort()
        signs = [sign(v) for _, v in values]
        flips = sum(
            a != b
            for a, b in zip(signs[:-1], signs[1:])
            if a != 0 and b != 0
        )
        trajectories.append(
            {
                "branch_index": int(key[0]),
                "requested_x": float(key[1]),
                "kind": key[2],
                "signs": [
                    {
                        "drive_scale": float(s),
                        "gradient": float(v),
                        "sign": sign(v),
                    }
                    for s, v in values
                ],
                "sign_flip_count": int(flips),
                "changes_sign": bool(flips > 0),
            }
        )

    def record_at(branch_index: int, scale: float) -> dict:
        return next(
            r
            for r in records
            if r["branch_index"] == branch_index
            and r["drive_scale"] == scale
        )

    baseline_proximal_length = []
    for bi in range(len(branches)):
        rec = record_at(bi, 1.0)
        row = next(
            x
            for x in rec["rows"]
            if x["requested_x"] == 0.10
            and x["kind"] == "length"
        )
        baseline_proximal_length.append(row)

    high_drive_diameter = []
    for bi in range(len(branches)):
        rec = record_at(bi, 3.0)
        high_drive_diameter.extend(
            [
                x
                for x in rec["rows"]
                if x["kind"] == "diameter"
            ]
        )

    low_drive_diameter = []
    for bi in range(len(branches)):
        rec = record_at(bi, 0.25)
        low_drive_diameter.extend(
            [
                x
                for x in rec["rows"]
                if x["kind"] == "diameter"
            ]
        )

    length_trajectories = [
        x for x in trajectories if x["kind"] == "length"
    ]
    diameter_trajectories = [
        x for x in trajectories if x["kind"] == "diameter"
    ]
    all_converged = all(r["all_converged"] for r in records)

    summary = {
        "object": (
            "within-cell branch panel for operating-point-dependent "
            "geometry-gradient signs"
        ),
        "fci_commit": FCI_COMMIT,
        "cell": {
            "morphology_identifier": model[
                "morphology_identifier"
            ],
            "species": model["species"],
            "layer": model["layer"],
            "compartments": int(len(graph["nodes"])),
            "branches": int(len(branches)),
        },
        "protocol": {
            "branches": [
                section_name(branch) for branch in branches
            ],
            "probe_x": list(PROBE_X),
            "nonlinear_site_x": list(SITE_X),
            "drive_scales": list(DRIVE_SCALES),
            "directions_per_branch_per_drive": int(
                2 * len(PROBE_X)
            ),
            "operating_points": int(
                len(branches) * len(DRIVE_SCALES)
            ),
        },
        "aggregate": {
            "gradient_trajectories": int(len(trajectories)),
            "trajectories_changing_sign": int(
                sum(x["changes_sign"] for x in trajectories)
            ),
            "length_trajectories_changing_sign": int(
                sum(
                    x["changes_sign"]
                    for x in length_trajectories
                )
            ),
            "diameter_trajectories_changing_sign": int(
                sum(
                    x["changes_sign"]
                    for x in diameter_trajectories
                )
            ),
            "baseline_1x_proximal_length_positive_branches": int(
                sum(x["sign"] > 0 for x in baseline_proximal_length)
            ),
            "baseline_1x_proximal_length_negative_branches": int(
                sum(x["sign"] < 0 for x in baseline_proximal_length)
            ),
            "low_025x_diameter_positive_directions": int(
                sum(x["sign"] > 0 for x in low_drive_diameter)
            ),
            "low_025x_diameter_negative_directions": int(
                sum(x["sign"] < 0 for x in low_drive_diameter)
            ),
            "high_3x_diameter_positive_directions": int(
                sum(x["sign"] > 0 for x in high_drive_diameter)
            ),
            "high_3x_diameter_negative_directions": int(
                sum(x["sign"] < 0 for x in high_drive_diameter)
            ),
            "duplicate_probe_nodes": duplicate_probe_nodes,
        },
        "trajectories": trajectories,
        "records": records,
        "all_converged": bool(all_converged),
        "classification": (
            "WITHIN_CELL_GRADIENT_REGIME_PANEL_VALID"
            if all_converged
            else "WITHIN_CELL_GRADIENT_REGIME_PANEL_FAILED"
        ),
        "stopping_line": (
            "This establishes or rejects within-cell branch recurrence "
            "of the drive-dependent sign regime. It is not a cross-cell "
            "generalization claim."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    a = summary["aggregate"]
    print("Operaattori within-cell gradient regime panel")
    print(
        f"cell {model['morphology_identifier']} "
        f"branches={len(branches)} "
        f"operating_points={summary['protocol']['operating_points']}"
    )
    print(
        "sign-changing trajectories: "
        f"{a['trajectories_changing_sign']}/"
        f"{a['gradient_trajectories']} "
        f"(length {a['length_trajectories_changing_sign']}/"
        f"{len(length_trajectories)}, diameter "
        f"{a['diameter_trajectories_changing_sign']}/"
        f"{len(diameter_trajectories)})"
    )
    print(
        "1x proximal length branches +/−: "
        f"{a['baseline_1x_proximal_length_positive_branches']}/"
        f"{a['baseline_1x_proximal_length_negative_branches']}"
    )
    print(
        "diameter directions at 0.25x +/−: "
        f"{a['low_025x_diameter_positive_directions']}/"
        f"{a['low_025x_diameter_negative_directions']}"
    )
    print(
        "diameter directions at 3x +/−: "
        f"{a['high_3x_diameter_positive_directions']}/"
        f"{a['high_3x_diameter_negative_directions']}"
    )
    print(f"classification: {summary['classification']}")
    if not all_converged:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
