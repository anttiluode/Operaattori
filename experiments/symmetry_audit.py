from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from gate10_real_neuron_scaffold import SOURCE_NAME, download_source
from gate22_embedding_vs_metric import (
    path_segments,
    relative_vector_difference,
    stretch_subtree_metric,
    transfer_signature,
)
from operaattori.cable_path import PassiveCableParams
from operaattori.real_scaffold import (
    build_matrix_scaffold,
    child_counts,
    descendant_counts,
    descendant_mask,
    edge_lengths,
    load_morphio_tree,
    reconstruct,
    twist_scaffold,
)


FREQUENCIES_HZ = np.asarray([1.0, 5.0, 20.0, 80.0, 300.0])
AXES = ("x", "y", "z")
ANGLES_DEG = (-35.0, -15.0, 15.0, 35.0)
METRIC_SCALE = 1.20


def affected_paths(tree, positions, pivot: int) -> list[dict]:
    mask = descendant_mask(tree.parents, pivot)
    cc = child_counts(tree.parents)
    rows: list[dict] = []

    for leaf in range(1, len(tree.parents)):
        if cc[leaf] != 0 or int(tree.section_types[leaf]) == 2 or not mask[leaf]:
            continue
        lengths, radii, nodes = path_segments(
            positions, tree.radii, tree.parents, leaf
        )
        if len(lengths) < 12:
            continue
        affected = np.asarray(
            [
                bool(mask[node] and mask[int(tree.parents[node])])
                for node in nodes
            ],
            dtype=bool,
        )
        affected_length = float(np.sum(lengths[affected]))
        if affected_length < 15.0:
            continue
        rows.append(
            {
                "leaf": int(leaf),
                "affected_length_um": affected_length,
                "affected_fraction": float(
                    affected_length / (float(np.sum(lengths)) + 1e-30)
                ),
            }
        )

    rows.sort(
        key=lambda r: (r["affected_length_um"], r["affected_fraction"]),
        reverse=True,
    )
    return rows


def choose_diverse_pivots(tree, positions, count: int, paths_per_pivot: int):
    cc = child_counts(tree.parents)
    dc = descendant_counts(tree.parents)
    n = len(tree.parents)
    viable = []

    for pivot in range(1, n):
        if (
            cc[pivot] < 2
            or dc[pivot] < 50
            or dc[pivot] > 0.40 * n
            or int(tree.section_types[pivot]) == 2
        ):
            continue
        paths = affected_paths(tree, positions, pivot)
        if len(paths) >= paths_per_pivot:
            viable.append((int(dc[pivot]), int(pivot), paths))

    if len(viable) < count:
        raise RuntimeError(
            f"only {len(viable)} viable non-axonal bifurcations for {count} requested"
        )

    viable.sort(key=lambda x: x[0])
    indices = np.rint(np.linspace(0, len(viable) - 1, count)).astype(int)
    chosen = []
    seen = set()
    for idx in indices:
        desc, pivot, paths = viable[int(idx)]
        if pivot in seen:
            continue
        seen.add(pivot)
        chosen.append((desc, pivot, paths[:paths_per_pivot]))

    if len(chosen) < count:
        for desc, pivot, paths in viable:
            if pivot not in seen:
                chosen.append((desc, pivot, paths[:paths_per_pivot]))
                seen.add(pivot)
            if len(chosen) == count:
                break

    return chosen[:count]


def path_transfer_difference(tree, base_pos, other_pos, leaf, params):
    l0, r0, _ = path_segments(base_pos, tree.radii, tree.parents, leaf)
    l1, r1, _ = path_segments(other_pos, tree.radii, tree.parents, leaf)
    if len(l0) != len(l1) or not np.allclose(r0, r1, rtol=0, atol=0):
        raise RuntimeError("path identity/radius changed during symmetry audit")

    z0, g0 = transfer_signature(l0, r0, FREQUENCIES_HZ, params)
    z1, g1 = transfer_signature(l1, r1, FREQUENCIES_HZ, params)
    dz = relative_vector_difference(z0, z1)
    dg = relative_vector_difference(g0, g1)
    return {
        "max": float(max(np.max(dz), np.max(dg))),
        "median": float(max(np.median(dz), np.median(dg))),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--morphology", type=Path)
    ap.add_argument("--pivots", type=int, default=6)
    ap.add_argument("--paths-per-pivot", type=int, default=3)
    ap.add_argument(
        "--out-dir", type=Path, default=Path("results/symmetry_audit")
    )
    args = ap.parse_args()

    if args.pivots < 3:
        raise ValueError("audit requires at least three pivots")
    if args.paths_per_pivot < 2:
        raise ValueError("audit requires at least two paths per pivot")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(out / "_source" / SOURCE_NAME)

    tree = load_morphio_tree(source)
    scaffold = build_matrix_scaffold(tree)
    base_pos, _ = reconstruct(scaffold)
    base_lengths = edge_lengths(base_pos, tree.parents)
    params = PassiveCableParams()

    chosen = choose_diverse_pivots(
        tree, base_pos, args.pivots, args.paths_per_pivot
    )

    bend_rows = []
    metric_rows = []

    for descendants, pivot, path_rows in chosen:
        mask = descendant_mask(tree.parents, pivot)
        leaves = [int(r["leaf"]) for r in path_rows]

        # The positive control changes only intrinsic translations inside the
        # same selected subtree.
        stretched = stretch_subtree_metric(scaffold, pivot, METRIC_SCALE)
        stretch_pos, _ = reconstruct(stretched)
        stretch_lengths = edge_lengths(stretch_pos, tree.parents)
        stretch_path_effects = [
            path_transfer_difference(
                tree, base_pos, stretch_pos, leaf, params
            )["median"]
            for leaf in leaves
        ]
        metric_rows.append(
            {
                "pivot": int(pivot),
                "descendants": int(descendants),
                "paths": leaves,
                "max_cable_length_change_um": float(
                    np.max(np.abs(stretch_lengths - base_lengths))
                ),
                "median_path_transfer_change": float(
                    np.median(stretch_path_effects)
                ),
                "min_path_transfer_change": float(
                    np.min(stretch_path_effects)
                ),
            }
        )

        for axis in AXES:
            for angle in ANGLES_DEG:
                bent = twist_scaffold(
                    scaffold,
                    pivot,
                    angle_degrees=float(angle),
                    axis=axis,
                )
                bent_pos, _ = reconstruct(bent)
                bent_lengths = edge_lengths(bent_pos, tree.parents)
                displacement = np.linalg.norm(bent_pos - base_pos, axis=1)
                path_effects = [
                    path_transfer_difference(
                        tree, base_pos, bent_pos, leaf, params
                    )["max"]
                    for leaf in leaves
                ]
                bend_rows.append(
                    {
                        "pivot": int(pivot),
                        "descendants": int(descendants),
                        "axis": axis,
                        "angle_degrees": float(angle),
                        "paths": leaves,
                        "max_distal_displacement_um": float(
                            np.max(displacement[mask])
                        ),
                        "max_cable_length_change_um": float(
                            np.max(np.abs(bent_lengths - base_lengths))
                        ),
                        "max_passive_transfer_change": float(
                            np.max(path_effects)
                        ),
                    }
                )

    iso_length = np.asarray(
        [r["max_cable_length_change_um"] for r in bend_rows], dtype=float
    )
    iso_transfer = np.asarray(
        [r["max_passive_transfer_change"] for r in bend_rows], dtype=float
    )
    iso_displacement = np.asarray(
        [r["max_distal_displacement_um"] for r in bend_rows], dtype=float
    )
    metric_effect = np.asarray(
        [r["median_path_transfer_change"] for r in metric_rows], dtype=float
    )

    aggregate = {
        "pivots": int(len(chosen)),
        "isometric_bend_trials": int(len(bend_rows)),
        "axes": list(AXES),
        "angles_degrees": list(ANGLES_DEG),
        "frequencies_hz": FREQUENCIES_HZ.tolist(),
        "paths_per_pivot": int(args.paths_per_pivot),
        "max_isometric_cable_length_change_um": float(np.max(iso_length)),
        "max_isometric_passive_transfer_change": float(np.max(iso_transfer)),
        "median_isometric_displacement_um": float(np.median(iso_displacement)),
        "max_isometric_displacement_um": float(np.max(iso_displacement)),
        "metric_stretch_factor": METRIC_SCALE,
        "median_metric_control_transfer_change": float(
            np.median(metric_effect)
        ),
        "fraction_metric_controls_over_1pct": float(
            np.mean(metric_effect >= 0.01)
        ),
    }

    passed = (
        len(chosen) >= args.pivots
        and len(bend_rows) == args.pivots * len(AXES) * len(ANGLES_DEG)
        and aggregate["max_isometric_cable_length_change_um"] < 1e-7
        and aggregate["max_isometric_passive_transfer_change"] < 1e-9
        and aggregate["max_isometric_displacement_um"] > 10.0
        and aggregate["median_metric_control_transfer_change"] >= 0.01
        and aggregate["fraction_metric_controls_over_1pct"] >= 0.80
    )

    classification = (
        "KNOWN_REEMBEDDING_SYMMETRY_REPLICATED_ON_REAL_SCAFFOLD"
        if passed
        else "SYMMETRY_AUDIT_FAILED"
    )

    summary = {
        "object": (
            "multi-pivot audit of classical cable invariance under local "
            "isometric re-embedding of real human cell 1125"
        ),
        "literature_fence": (
            "Lopez-Sanchez & Romero, Phys Rev E 95, 022403 (2017): for a "
            "constant circular cross-section their generalized cable equation "
            "depends on neither curvature nor torsion. This audit is therefore "
            "a software/causal replication, not a novelty claim."
        ),
        "protocol": {
            "pivot_selection": (
                "deterministic quantiles across viable non-axonal bifurcations "
                "with >=50 descendants and >= requested affected tip paths"
            ),
            "isometric_axes": list(AXES),
            "isometric_angles_degrees": list(ANGLES_DEG),
            "metric_positive_control": "20% intrinsic subtree stretch",
            "thresholds_locked_before_run": {
                "isometric_cable_length_change_um_max": 1e-7,
                "isometric_passive_transfer_change_max": 1e-9,
                "max_visual_displacement_um_min": 10.0,
                "metric_median_transfer_change_min": 0.01,
                "metric_controls_over_1pct_fraction_min": 0.80,
            },
        },
        "aggregate": aggregate,
        "classification": classification,
        "pivots": [
            {
                "pivot": int(pivot),
                "descendants": int(descendants),
                "paths": paths,
            }
            for descendants, pivot, paths in chosen
        ],
        "isometric_trials": bend_rows,
        "metric_controls": metric_rows,
    }

    (out / "symmetry_audit.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print("Operaattori symmetry audit — intrinsic cable vs embedding")
    print()
    print(f"pivots:                              {aggregate['pivots']}")
    print(
        "isometric bend trials:               "
        f"{aggregate['isometric_bend_trials']}"
    )
    print(
        "max visual displacement:              "
        f"{aggregate['max_isometric_displacement_um']:.3f} um"
    )
    print(
        "max isometric cable-length change:    "
        f"{aggregate['max_isometric_cable_length_change_um']:.3e} um"
    )
    print(
        "max isometric passive-transfer change:"
        f" {aggregate['max_isometric_passive_transfer_change']:.3e}"
    )
    print(
        "median 20% metric-control change:      "
        f"{aggregate['median_metric_control_transfer_change']:.4f}"
    )
    print(
        "metric controls >1%:                   "
        f"{aggregate['fraction_metric_controls_over_1pct']:.3f}"
    )
    print()
    print(f"classification: {classification}")

    assert passed, classification


if __name__ == "__main__":
    main()
