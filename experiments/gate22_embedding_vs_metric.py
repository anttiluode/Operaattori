from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from gate10_real_neuron_scaffold import SOURCE_NAME, download_source
from operaattori.cable_path import (
    PassiveCableParams,
    path_abcd,
    relative_complex_difference,
    sealed_distal_voltage_gain,
    sealed_input_impedance,
)
from operaattori.real_scaffold import (
    MatrixScaffold,
    build_matrix_scaffold,
    child_counts,
    choose_twist_pivot,
    descendant_mask,
    edge_lengths,
    load_morphio_tree,
    reconstruct,
    twist_scaffold,
)


def path_nodes_to_root(parents: np.ndarray, leaf: int) -> list[int]:
    out = [int(leaf)]
    i = int(leaf)
    while i > 0:
        i = int(parents[i])
        out.append(i)
    out.reverse()
    return out


def path_segments(
    positions: np.ndarray,
    radii_nodes: np.ndarray,
    parents: np.ndarray,
    leaf: int,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    nodes = path_nodes_to_root(parents, leaf)
    lengths: list[float] = []
    radii: list[float] = []
    child_nodes: list[int] = []
    for child in nodes[1:]:
        parent = int(parents[child])
        # Synthetic soma-centroid -> first-neurite edge is scaffold
        # bookkeeping rather than a cylindrical cable segment.
        if parent == 0:
            continue
        length = float(
            np.linalg.norm(positions[child] - positions[parent])
        )
        if length <= 1e-8:
            continue
        radius = 0.5 * float(
            radii_nodes[child] + radii_nodes[parent]
        )
        radius = max(radius, 0.15)
        lengths.append(length)
        radii.append(radius)
        child_nodes.append(int(child))
    return (
        np.asarray(lengths, dtype=float),
        np.asarray(radii, dtype=float),
        child_nodes,
    )


def stretch_subtree_metric(
    scaffold: MatrixScaffold,
    pivot: int,
    scale: float,
) -> MatrixScaffold:
    if scale <= 0:
        raise ValueError("scale must be positive")
    mask = descendant_mask(scaffold.parents, pivot)
    local = scaffold.local_transforms.copy()

    # Keep the parent->pivot attachment untouched. Scale only edges whose two
    # endpoints are inside the selected subtree.
    for i in range(1, len(scaffold.parents)):
        p = int(scaffold.parents[i])
        if mask[i] and mask[p]:
            local[i, :3, 3] *= float(scale)

    out = replace(scaffold, local_transforms=local)
    out.validate()
    return out


def transfer_signature(
    lengths: np.ndarray,
    radii: np.ndarray,
    frequencies: np.ndarray,
    params: PassiveCableParams,
) -> tuple[np.ndarray, np.ndarray]:
    z = []
    gain = []
    for f in frequencies:
        M = path_abcd(lengths, radii, float(f), params)
        z.append(sealed_input_impedance(M))
        gain.append(sealed_distal_voltage_gain(M))
    return (
        np.asarray(z, dtype=np.complex128),
        np.asarray(gain, dtype=np.complex128),
    )


def relative_vector_difference(
    a: np.ndarray,
    b: np.ndarray,
) -> np.ndarray:
    return np.asarray(
        [
            relative_complex_difference(complex(x), complex(y))
            for x, y in zip(a, b)
        ],
        dtype=float,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--morphology", type=Path)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/gate22"),
    )
    ap.add_argument("--paths", type=int, default=32)
    ap.add_argument("--twist-deg", type=float, default=35.0)
    ap.add_argument("--stretch", type=float, default=1.20)
    args = ap.parse_args()

    if abs(args.twist_deg - 35.0) > 1e-12:
        raise ValueError("Gate 22 twist is locked to 35 degrees")
    if abs(args.stretch - 1.20) > 1e-12:
        raise ValueError("Gate 22 metric stretch is locked to 1.20")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(
        out / "_source" / SOURCE_NAME
    )

    tree = load_morphio_tree(source)
    scaffold = build_matrix_scaffold(tree)
    original_pos, _ = reconstruct(scaffold)

    pivot = choose_twist_pivot(tree.parents)
    mask = descendant_mask(tree.parents, pivot)

    twisted = twist_scaffold(
        scaffold,
        pivot,
        angle_degrees=args.twist_deg,
        axis="y",
    )
    twisted_pos, _ = reconstruct(twisted)

    stretched = stretch_subtree_metric(
        scaffold,
        pivot,
        args.stretch,
    )
    stretched_pos, _ = reconstruct(stretched)

    original_lengths_all = edge_lengths(original_pos, tree.parents)
    twist_lengths_all = edge_lengths(twisted_pos, tree.parents)
    stretch_lengths_all = edge_lengths(stretched_pos, tree.parents)

    displacement = np.linalg.norm(twisted_pos - original_pos, axis=1)

    cc = child_counts(tree.parents)
    leaves = [
        i
        for i in range(1, len(tree.parents))
        if cc[i] == 0
        and int(tree.section_types[i]) != 2
        and mask[i]
    ]

    candidates = []
    for leaf in leaves:
        l0, r0, nodes = path_segments(
            original_pos,
            tree.radii,
            tree.parents,
            leaf,
        )
        if len(l0) < 20:
            continue
        affected = np.asarray(
            [
                mask[node] and mask[int(tree.parents[node])]
                for node in nodes
            ],
            dtype=bool,
        )
        affected_length = float(np.sum(l0[affected]))
        total_length = float(np.sum(l0))
        frac = affected_length / (total_length + 1e-30)
        if affected_length < 20.0:
            continue
        candidates.append(
            (frac, affected_length, total_length, int(leaf))
        )

    candidates.sort(reverse=True)
    selected = candidates[: min(args.paths, len(candidates))]
    if len(selected) < 8:
        raise RuntimeError("too few affected dendritic tip paths")

    frequencies = np.geomspace(1.0, 300.0, 36)
    params = PassiveCableParams()
    rows = []

    for frac, affected_length, total_length, leaf in selected:
        l0, r0, _ = path_segments(
            original_pos, tree.radii, tree.parents, leaf
        )
        lt, rt, _ = path_segments(
            twisted_pos, tree.radii, tree.parents, leaf
        )
        ls, rs, _ = path_segments(
            stretched_pos, tree.radii, tree.parents, leaf
        )

        if not (
            len(l0) == len(lt)
            and len(l0) == len(ls)
            and np.allclose(r0, rt, rtol=0, atol=0)
            and np.allclose(r0, rs, rtol=0, atol=0)
        ):
            raise RuntimeError("Gate-22 path identity changed")

        z0, g0 = transfer_signature(l0, r0, frequencies, params)
        zt, gt = transfer_signature(lt, rt, frequencies, params)
        zs, gs = transfer_signature(ls, rs, frequencies, params)

        twist_z = relative_vector_difference(z0, zt)
        twist_g = relative_vector_difference(g0, gt)
        stretch_z = relative_vector_difference(z0, zs)
        stretch_g = relative_vector_difference(g0, gs)

        rows.append(
            {
                "leaf": int(leaf),
                "segments": int(len(l0)),
                "total_path_length_um": total_length,
                "affected_length_um": affected_length,
                "affected_fraction": frac,
                "twist": {
                    "max_impedance_relative_difference": float(
                        np.max(twist_z)
                    ),
                    "max_gain_relative_difference": float(
                        np.max(twist_g)
                    ),
                    "median_impedance_relative_difference": float(
                        np.median(twist_z)
                    ),
                    "median_gain_relative_difference": float(
                        np.median(twist_g)
                    ),
                },
                "metric_stretch": {
                    "median_impedance_relative_difference": float(
                        np.median(stretch_z)
                    ),
                    "median_gain_relative_difference": float(
                        np.median(stretch_g)
                    ),
                    "max_impedance_relative_difference": float(
                        np.max(stretch_z)
                    ),
                    "max_gain_relative_difference": float(
                        np.max(stretch_g)
                    ),
                },
            }
        )

    twist_max = max(
        max(
            r["twist"]["max_impedance_relative_difference"],
            r["twist"]["max_gain_relative_difference"],
        )
        for r in rows
    )
    stretch_path = np.asarray(
        [
            max(
                r["metric_stretch"][
                    "median_impedance_relative_difference"
                ],
                r["metric_stretch"]["median_gain_relative_difference"],
            )
            for r in rows
        ],
        dtype=float,
    )

    max_twist_length_change = float(
        np.max(np.abs(twist_lengths_all - original_lengths_all))
    )

    aggregate = {
        "paths": int(len(rows)),
        "pivot_node": int(pivot),
        "twist_degrees": float(args.twist_deg),
        "stretch_factor": float(args.stretch),
        "max_distal_3d_displacement_um": float(
            np.max(displacement[mask])
        ),
        "max_outside_3d_displacement_um": float(
            np.max(displacement[~mask])
        ),
        "max_twist_cable_length_change_um": max_twist_length_change,
        "max_twist_passive_transfer_relative_difference": float(
            twist_max
        ),
        "median_metric_stretch_transfer_relative_difference": float(
            np.median(stretch_path)
        ),
        "fraction_metric_stretch_paths_over_1pct": float(
            np.mean(stretch_path >= 0.01)
        ),
        "total_cable_length_original_um": float(
            np.sum(original_lengths_all)
        ),
        "total_cable_length_twisted_um": float(
            np.sum(twist_lengths_all)
        ),
        "total_cable_length_stretched_um": float(
            np.sum(stretch_lengths_all)
        ),
    }

    if (
        aggregate["max_distal_3d_displacement_um"] > 10.0
        and max_twist_length_change < 1e-7
        and twist_max < 1e-9
        and aggregate[
            "median_metric_stretch_transfer_relative_difference"
        ] >= 0.01
    ):
        classification = (
            "CABLE_MODEL_IGNORES_ISOMETRIC_3D_EMBEDDING"
        )
        interpretation = (
            "A large parent-local SE(3) bend moves the real distal arbor "
            "without changing its cable metric and leaves the classical "
            "passive electrical operator invariant to numerical precision. "
            "Changing the intrinsic metric on the same subtree does change "
            "transport. XYZ embedding is therefore not itself an internal "
            "cable-computation degree of freedom in this model."
        )
    elif twist_max >= 1e-9:
        classification = "UNEXPECTED_EMBEDDING_DEPENDENCE"
        interpretation = (
            "The supposedly isometric scaffold bend changed passive cable "
            "transfer beyond the preregistered numerical tolerance."
        )
    else:
        classification = "METRIC_POSITIVE_CONTROL_TOO_WEAK"
        interpretation = (
            "The isometric bend is electrically invariant as expected, but "
            "the locked 20% intrinsic metric stretch did not clear the 1% "
            "positive-control threshold."
        )

    summary = {
        "gate": 22,
        "object": (
            "extrinsic XYZ embedding versus intrinsic cable metric on the "
            "real cell-1125 matrix scaffold"
        ),
        "protocol": {
            "same_gate10_pivot_rule": True,
            "twist_degrees": float(args.twist_deg),
            "metric_stretch_factor": float(args.stretch),
            "radii_unchanged": True,
            "topology_unchanged": True,
            "frequencies_hz": frequencies.tolist(),
            "dendritic_paths": int(len(rows)),
            "electrical_model": (
                "same classical passive ABCD cable matrices as Gate 11"
            ),
            "thresholds_locked_before_run": {
                "distal_displacement_um_min": 10.0,
                "twist_length_change_um_max": 1e-7,
                "twist_transfer_relative_difference_max": 1e-9,
                "metric_stretch_median_transfer_difference_min": 0.01,
            },
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "paths": rows,
        "stopping_line": (
            "If the isometric-embedding invariance passes, 3-D scaffold "
            "geometry must couple to something spatial outside the intrinsic "
            "cable equation before a pure bend can compute differently. The "
            "next honest test is a fixed structured external field sampled "
            "by the bent versus unbent real arbor."
        ),
    }

    (out / "gate22_embedding_vs_metric.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Operaattori Gate 22 — scaffold embedding versus cable metric")
    print()
    print(f"paths:                              {aggregate['paths']}")
    print(f"pivot:                              {pivot}")
    print(
        "max distal 3-D displacement:       "
        f"{aggregate['max_distal_3d_displacement_um']:.3f} um"
    )
    print(
        "max twist cable-length change:     "
        f"{max_twist_length_change:.3e} um"
    )
    print(
        "max twist transfer difference:     "
        f"{twist_max:.3e}"
    )
    print(
        "median metric-stretch difference:  "
        f"{aggregate['median_metric_stretch_transfer_relative_difference']:.4f}"
    )
    print(
        "metric-stretch paths >1%:          "
        f"{aggregate['fraction_metric_stretch_paths_over_1pct']:.3f}"
    )
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(rows) >= 8
    assert np.all(np.isfinite(stretch_path))
    assert classification in {
        "CABLE_MODEL_IGNORES_ISOMETRIC_3D_EMBEDDING",
        "UNEXPECTED_EMBEDDING_DEPENDENCE",
        "METRIC_POSITIVE_CONTROL_TOO_WEAK",
    }


if __name__ == "__main__":
    main()
