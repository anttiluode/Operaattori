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
from gradient_operating_point import conductance_program
from gradient_transport_feedback import (
    cosine,
    sign,
    simulate_decomposition,
)
from real_metric_tangent import metric_tangent

DRIVE_SCALES = (0.5, 1.0, 2.0)


def site_directions(graph: dict, branch) -> tuple[list, list[dict]]:
    directions = []
    meta = []
    for site_index, x in enumerate(SITE_X):
        node_index = int(
            node_for_section_x(graph, branch, float(x))
        )
        node = graph["nodes"][node_index]
        for kind in ("length", "diameter"):
            directions.append(
                metric_tangent(graph, node_index, kind)
            )
            meta.append(
                {
                    "site_index": int(site_index),
                    "site_x": float(x),
                    "node": node_index,
                    "kind": kind,
                    "actual_x": float(node["x"]),
                }
            )
    return directions, meta


def summarize_scale(records: list[dict], scale: float) -> dict:
    rows = [x for x in records if x["drive_scale"] == scale]
    ratios = np.asarray(
        [x["feedback_to_transport_norm_ratio"] for x in rows],
        dtype=float,
    )
    cosines = np.asarray(
        [x["full_vs_transport_cosine"] for x in rows],
        dtype=float,
    )
    local = np.asarray(
        [x["max_local_depol_mV"] for x in rows],
        dtype=float,
    )
    return {
        "drive_scale": float(scale),
        "cells": int(len(rows)),
        "median_max_local_depol_mV": float(np.median(local)),
        "median_full_vs_transport_cosine": float(
            np.median(cosines)
        ),
        "min_full_vs_transport_cosine": float(np.min(cosines)),
        "median_feedback_to_transport_norm_ratio": float(
            np.median(ratios)
        ),
        "cells_feedback_norm_gt_transport": int(
            np.sum(ratios > 1.0)
        ),
        "total_sign_flips_full_vs_transport": int(
            sum(x["sign_flips"] for x in rows)
        ),
        "length_sign_flips_full_vs_transport": int(
            sum(x["length_sign_flips"] for x in rows)
        ),
        "diameter_sign_flips_full_vs_transport": int(
            sum(x["diameter_sign_flips"] for x in rows)
        ),
        "full_diameter_positive": int(
            sum(x["full_diameter_positive"] for x in rows)
        ),
        "transport_diameter_positive": int(
            sum(x["transport_diameter_positive"] for x in rows)
        ),
        "full_length_positive": int(
            sum(x["full_length_positive"] for x in rows)
        ),
        "transport_length_positive": int(
            sum(x["transport_length_positive"] for x in rows)
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/"
            "cross_cell_gradient_decomposition.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source not pinned")

    setup_neuron(fci_root)
    panel = load_panel(fci_root)
    records = []

    for ci, model in enumerate(panel):
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
        soma_node = node_for_section_x(
            graph, cell.soma[0], 0.5
        )
        directions, meta = site_directions(graph, branch)

        for scale in DRIVE_SCALES:
            ga, gn = conductance_program(scale)
            out = simulate_decomposition(
                graph["G_uS"],
                graph["C_nF"],
                site_nodes,
                soma_node,
                ga,
                gn,
                directions,
            )
            peak_index = int(np.argmax(out["soma_mV"]))
            full = out["full_tangent"][:, peak_index]
            transport = out[
                "frozen_current_tangent"
            ][:, peak_index]
            feedback = out[
                "feedback_tangent"
            ][:, peak_index]

            direction_rows = []
            for mi, gf, gt, gb in zip(
                meta, full, transport, feedback
            ):
                direction_rows.append(
                    {
                        **mi,
                        "full": float(gf),
                        "transport": float(gt),
                        "feedback": float(gb),
                        "full_sign": sign(gf),
                        "transport_sign": sign(gt),
                        "sign_flip": bool(
                            sign(gf) != 0
                            and sign(gt) != 0
                            and sign(gf) != sign(gt)
                        ),
                    }
                )

            length_rows = [
                x for x in direction_rows
                if x["kind"] == "length"
            ]
            diameter_rows = [
                x for x in direction_rows
                if x["kind"] == "diameter"
            ]

            transport_norm = float(
                np.linalg.norm(transport)
            )
            feedback_norm = float(
                np.linalg.norm(feedback)
            )
            records.append(
                {
                    "cell_order": int(model["order"]),
                    "species": model["species"],
                    "layer": model["layer"],
                    "morphology_identifier": model[
                        "morphology_identifier"
                    ],
                    "compartments": int(len(graph["nodes"])),
                    "branch": section_name(branch),
                    "drive_scale": float(scale),
                    "peak_mV": float(
                        out["soma_mV"][peak_index]
                    ),
                    "max_local_depol_mV": float(
                        np.max(out["local_depol_mV"])
                    ),
                    "full_l2": float(
                        np.linalg.norm(full)
                    ),
                    "transport_l2": transport_norm,
                    "feedback_l2": feedback_norm,
                    "feedback_to_transport_norm_ratio": float(
                        feedback_norm
                        / (transport_norm + 1e-30)
                    ),
                    "full_vs_transport_cosine": cosine(
                        full, transport
                    ),
                    "sign_flips": int(
                        sum(x["sign_flip"] for x in direction_rows)
                    ),
                    "length_sign_flips": int(
                        sum(x["sign_flip"] for x in length_rows)
                    ),
                    "diameter_sign_flips": int(
                        sum(x["sign_flip"] for x in diameter_rows)
                    ),
                    "full_length_positive": int(
                        sum(x["full_sign"] > 0 for x in length_rows)
                    ),
                    "transport_length_positive": int(
                        sum(
                            x["transport_sign"] > 0
                            for x in length_rows
                        )
                    ),
                    "full_diameter_positive": int(
                        sum(
                            x["full_sign"] > 0
                            for x in diameter_rows
                        )
                    ),
                    "transport_diameter_positive": int(
                        sum(
                            x["transport_sign"] > 0
                            for x in diameter_rows
                        )
                    ),
                    "directions": direction_rows,
                    "all_converged": bool(
                        out["all_converged"]
                    ),
                }
            )

        print(
            f"[{ci+1:02d}/24] {model['species']:5s} "
            f"{model['morphology_identifier']:>12s} "
            f"N={len(graph['nodes'])} "
            f"branch={section_name(branch)}"
        )

    all_converged = all(
        x["all_converged"] for x in records
    )
    by_scale = [
        summarize_scale(records, scale)
        for scale in DRIVE_SCALES
    ]

    ratios_by_cell = {}
    for model in panel:
        mid = model["morphology_identifier"]
        cell_rows = [
            x for x in records
            if x["morphology_identifier"] == mid
        ]
        ratios_by_cell[mid] = {
            "species": model["species"],
            "layer": model["layer"],
            "feedback_gt_transport_any_drive": bool(
                any(
                    x["feedback_to_transport_norm_ratio"] > 1.0
                    for x in cell_rows
                )
            ),
            "sign_flip_any_drive": bool(
                any(x["sign_flips"] > 0 for x in cell_rows)
            ),
            "records": [
                {
                    "drive_scale": x["drive_scale"],
                    "ratio": x[
                        "feedback_to_transport_norm_ratio"
                    ],
                    "cosine": x[
                        "full_vs_transport_cosine"
                    ],
                    "sign_flips": x["sign_flips"],
                }
                for x in cell_rows
            ],
        }

    summary = {
        "object": (
            "cross-cell recurrence of transport versus nonlinear-feedback "
            "geometry-gradient decomposition"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "cells": int(len(panel)),
            "branch_per_cell": 1,
            "sites_per_branch": len(SITE_X),
            "directions_per_cell_per_drive": 6,
            "drive_scales": list(DRIVE_SCALES),
            "operating_points": int(
                len(panel) * len(DRIVE_SCALES)
            ),
        },
        "aggregate_by_drive": by_scale,
        "cells_feedback_gt_transport_any_drive": int(
            sum(
                x["feedback_gt_transport_any_drive"]
                for x in ratios_by_cell.values()
            )
        ),
        "cells_with_any_full_transport_sign_flip": int(
            sum(
                x["sign_flip_any_drive"]
                for x in ratios_by_cell.values()
            )
        ),
        "per_cell": ratios_by_cell,
        "records": records,
        "classification": (
            "CROSS_CELL_TRANSPORT_NMDA_REGIME_PANEL_VALID"
            if all_converged
            else "CROSS_CELL_TRANSPORT_NMDA_REGIME_PANEL_FAILED"
        ),
        "stopping_line": (
            "This is a cross-cell sensitivity panel on one deterministic "
            "branch per morphology. It does not prove that a scalar local "
            "state predicts the decomposition or that the same pattern "
            "holds for every branch."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori cross-cell transport/NMDA gradient panel")
    for x in by_scale:
        print(
            f"drive {x['drive_scale']:.1f}: "
            f"median local={x['median_max_local_depol_mV']:.3g} mV "
            f"median cos={x['median_full_vs_transport_cosine']:+.3f} "
            f"median fb/trans={x['median_feedback_to_transport_norm_ratio']:.3f} "
            f"fb>trans cells={x['cells_feedback_norm_gt_transport']}/24 "
            f"signflips={x['total_sign_flips_full_vs_transport']}/144"
        )
    print(
        "cells feedback>transport at any drive: "
        f"{summary['cells_feedback_gt_transport_any_drive']}/24"
    )
    print(
        "cells with any full/transport sign flip: "
        f"{summary['cells_with_any_full_transport_sign_flip']}/24"
    )
    print(f"classification: {summary['classification']}")
    if not all_converged:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
