from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from gate10_real_neuron_scaffold import SOURCE_NAME, download_source
from gate13_branching_transport import build_electrical_tree
from operaattori.cable_path import PassiveCableParams
from operaattori.real_scaffold import (
    build_matrix_scaffold,
    choose_twist_pivot,
    descendant_mask,
    edge_lengths,
    load_morphio_tree,
    reconstruct,
    twist_scaffold,
)
from operaattori.tree_cable import solve_tree_frequency


def fibonacci_sphere(n: int) -> np.ndarray:
    i = np.arange(n, dtype=float)
    golden = np.pi * (3.0 - np.sqrt(5.0))
    y = 1.0 - 2.0 * (i + 0.5) / float(n)
    r = np.sqrt(np.maximum(0.0, 1.0 - y * y))
    theta = golden * i
    x = np.cos(theta) * r
    z = np.sin(theta) * r
    return np.column_stack([x, y, z])


def rel_complex(a: complex, b: complex) -> float:
    return float(
        abs(a - b) / (0.5 * (abs(a) + abs(b)) + 1e-30)
    )


def normalized_injection(
    positions: np.ndarray,
    pivot_position: np.ndarray,
    direction: np.ndarray,
    beta: float,
    spatial_scale: float,
    node_area: np.ndarray,
    sample_mask: np.ndarray,
) -> np.ndarray:
    proj = (
        (positions - pivot_position[None, :])
        @ np.asarray(direction, dtype=float)
    ) / float(spatial_scale)
    expo = np.clip(float(beta) * proj, -12.0, 12.0)
    field = np.exp(expo)
    raw = np.zeros(len(positions), dtype=float)
    raw[sample_mask] = field[sample_mask] * node_area[sample_mask]
    total = float(np.sum(raw))
    if total <= 0 or not np.isfinite(total):
        raise FloatingPointError("invalid Gate-23 spatial drive")
    # One arbitrary unit of total injected current for every condition.
    return raw / total


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--morphology", type=Path)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/gate23"),
    )
    ap.add_argument("--directions", type=int, default=48)
    ap.add_argument("--twist-deg", type=float, default=35.0)
    args = ap.parse_args()

    if args.directions != 48:
        raise ValueError("Gate 23 is locked to 48 field directions")
    if abs(args.twist_deg - 35.0) > 1e-12:
        raise ValueError("Gate 23 is locked to the Gate-22 35-degree bend")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(
        out / "_source" / SOURCE_NAME
    )
    tree = load_morphio_tree(source)

    scaffold = build_matrix_scaffold(tree)
    original_pos, _ = reconstruct(scaffold)
    pivot = choose_twist_pivot(tree.parents)
    material_mask = descendant_mask(tree.parents, pivot)

    bent = twist_scaffold(
        scaffold,
        pivot,
        angle_degrees=args.twist_deg,
        axis="y",
    )
    bent_pos, _ = reconstruct(bent)

    original_lengths_geom = edge_lengths(original_pos, tree.parents)
    bent_lengths_geom = edge_lengths(bent_pos, tree.parents)
    intrinsic_length_error = float(
        np.max(np.abs(original_lengths_geom - bent_lengths_geom))
    )

    parents, lengths, radii, active, clamped = build_electrical_tree(tree)
    sample_mask = material_mask & active & (~clamped)
    if int(np.sum(sample_mask)) < 100:
        raise RuntimeError("too few active material nodes in Gate-23 subtree")

    # Surface proxy for distributing a spatially sampled drive over physical
    # cable. It is identical before/after the isometric bend.
    node_area = np.zeros(len(tree.parents), dtype=float)
    node_area[sample_mask] = (
        2.0
        * np.pi
        * radii[sample_mask]
        * lengths[sample_mask]
    )

    pivot_position = original_pos[pivot].copy()
    rel = original_pos[sample_mask] - pivot_position[None, :]
    spatial_scale = float(
        np.sqrt(np.mean(np.sum(rel * rel, axis=1)))
    )
    if spatial_scale <= 1e-9:
        raise RuntimeError("degenerate Gate-23 spatial scale")

    directions = fibonacci_sphere(args.directions)
    betas = (0.5, 1.0, 2.0)
    frequencies = (5.0, 20.0, 80.0)
    params = PassiveCableParams()

    # The electrical transfer is computed once from the original intrinsic
    # cable. Gate 22 established that the isometric bend does not change it.
    transfer_by_freq: dict[float, np.ndarray] = {}
    for freq in frequencies:
        state = solve_tree_frequency(
            parents,
            lengths,
            radii,
            active,
            clamped,
            freq,
            params,
        )
        transfer_by_freq[freq] = state.transfer_to_clamp.copy()

    rows = []
    fixed_diffs = []
    material_locked_diffs = []

    for beta in betas:
        for di, direction in enumerate(directions):
            inj_original = normalized_injection(
                original_pos,
                pivot_position,
                direction,
                beta,
                spatial_scale,
                node_area,
                sample_mask,
            )
            inj_bent = normalized_injection(
                bent_pos,
                pivot_position,
                direction,
                beta,
                spatial_scale,
                node_area,
                sample_mask,
            )

            l1 = float(np.sum(inj_original))
            l2 = float(np.sum(inj_bent))
            if abs(l1 - 1.0) > 1e-12 or abs(l2 - 1.0) > 1e-12:
                raise RuntimeError("Gate-23 equal-total normalization failed")

            field_tv = float(
                0.5 * np.sum(np.abs(inj_original - inj_bent))
            )

            for freq in frequencies:
                T = transfer_by_freq[freq]
                y_original = complex(np.sum(T * inj_original))
                y_bent_fixed_world = complex(np.sum(T * inj_bent))

                # Material-locked attacker: attach the original sampled values
                # to the same material nodes after bending. The vector is
                # therefore exactly inj_original and the response should be
                # invariant.
                y_bent_material_locked = complex(
                    np.sum(T * inj_original)
                )

                fixed_diff = rel_complex(
                    y_original,
                    y_bent_fixed_world,
                )
                locked_diff = rel_complex(
                    y_original,
                    y_bent_material_locked,
                )
                fixed_diffs.append(fixed_diff)
                material_locked_diffs.append(locked_diff)

                rows.append(
                    {
                        "beta": float(beta),
                        "direction_index": int(di),
                        "frequency_hz": float(freq),
                        "equal_total_injection_original": l1,
                        "equal_total_injection_bent": l2,
                        "injection_total_variation": field_tv,
                        "fixed_world_output_relative_difference": fixed_diff,
                        "material_locked_output_relative_difference": locked_diff,
                        "original_output_real": float(np.real(y_original)),
                        "original_output_imag": float(np.imag(y_original)),
                        "bent_output_real": float(
                            np.real(y_bent_fixed_world)
                        ),
                        "bent_output_imag": float(
                            np.imag(y_bent_fixed_world)
                        ),
                    }
                )

    fixed = np.asarray(fixed_diffs, dtype=float)
    locked = np.asarray(material_locked_diffs, dtype=float)

    aggregate = {
        "sample_nodes": int(np.sum(sample_mask)),
        "field_directions": int(args.directions),
        "gradient_strengths": list(betas),
        "frequencies_hz": list(frequencies),
        "conditions": int(len(rows)),
        "spatial_scale_um": spatial_scale,
        "max_intrinsic_cable_length_change_um": intrinsic_length_error,
        "max_material_locked_output_relative_difference": float(
            np.max(locked)
        ),
        "median_fixed_world_output_relative_difference": float(
            np.median(fixed)
        ),
        "mean_fixed_world_output_relative_difference": float(
            np.mean(fixed)
        ),
        "fraction_fixed_world_conditions_over_5pct": float(
            np.mean(fixed >= 0.05)
        ),
        "fraction_fixed_world_conditions_over_10pct": float(
            np.mean(fixed >= 0.10)
        ),
        "max_fixed_world_output_relative_difference": float(
            np.max(fixed)
        ),
    }

    if (
        intrinsic_length_error < 1e-7
        and aggregate[
            "max_material_locked_output_relative_difference"
        ] < 1e-12
        and aggregate[
            "median_fixed_world_output_relative_difference"
        ] >= 0.05
        and aggregate[
            "fraction_fixed_world_conditions_over_5pct"
        ] >= 0.50
    ):
        classification = (
            "SPATIAL_COUPLING_MAKES_EMBEDDING_FUNCTIONAL"
        )
        interpretation = (
            "The same isometric matrix bend that is invisible to the "
            "intrinsic cable operator changes the equal-total somatic readout "
            "when material nodes sample a fixed structured world field. The "
            "material-locked control remains invariant, locating causality in "
            "the relative geometry between scaffold and environment."
        )
    elif aggregate[
        "max_material_locked_output_relative_difference"
    ] >= 1e-12:
        classification = "MATERIAL_LOCKED_INVARIANCE_FAILED"
        interpretation = (
            "The control that keeps field values attached to material nodes "
            "changed unexpectedly, so the spatial-coupling result is invalid."
        )
    else:
        classification = "FIXED_WORLD_FIELD_EFFECT_WEAK"
        interpretation = (
            "The isometric bend changes world-space sampling but the locked "
            "gradient family does not produce a robust >=5% somatic effect "
            "after equal-total normalization."
        )

    summary = {
        "gate": 23,
        "object": (
            "fixed world-space field coupled to an isometrically movable "
            "real-neuron scaffold"
        ),
        "protocol": {
            "same_gate22_isometric_bend": True,
            "twist_degrees": float(args.twist_deg),
            "intrinsic_tree_reused_for_both_embeddings": True,
            "field": (
                "positive exponential gradient exp(beta*d dot "
                "(x-xpivot)/R)"
            ),
            "directions": int(args.directions),
            "betas": list(betas),
            "frequencies_hz": list(frequencies),
            "surface_weighted_sampling": True,
            "equal_total_injection_per_embedding": True,
            "material_locked_control": True,
            "thresholds_locked_before_run": {
                "intrinsic_length_change_um_max": 1e-7,
                "material_locked_output_difference_max": 1e-12,
                "median_fixed_world_output_difference_min": 0.05,
                "fraction_fixed_world_conditions_over_5pct_min": 0.50,
            },
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "conditions": rows,
        "stopping_line": (
            "A positive spatial-coupling result earns an adaptive scaffold "
            "test: a local rule may modify geometry only if improvement is "
            "measured on held-out world fields and compared with simpler "
            "weight-only/input-gain attackers. It still does not earn a "
            "biological growth claim."
        ),
    }

    (out / "gate23_spatial_field.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Operaattori Gate 23 — fixed world field on movable scaffold")
    print()
    print(f"sample nodes:                         {aggregate['sample_nodes']}")
    print(f"field/frequency conditions:           {aggregate['conditions']}")
    print(
        "max intrinsic length change:         "
        f"{intrinsic_length_error:.3e} um"
    )
    print(
        "material-locked max output diff:      "
        f"{aggregate['max_material_locked_output_relative_difference']:.3e}"
    )
    print(
        "fixed-world median output diff:       "
        f"{aggregate['median_fixed_world_output_relative_difference']:.4f}"
    )
    print(
        "fixed-world conditions >5%:           "
        f"{aggregate['fraction_fixed_world_conditions_over_5pct']:.3f}"
    )
    print(
        "fixed-world conditions >10%:          "
        f"{aggregate['fraction_fixed_world_conditions_over_10pct']:.3f}"
    )
    print(
        "fixed-world max output diff:          "
        f"{aggregate['max_fixed_world_output_relative_difference']:.4f}"
    )
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(rows) == 48 * 3 * 3
    assert np.all(np.isfinite(fixed))
    assert np.all(np.isfinite(locked))
    assert classification in {
        "SPATIAL_COUPLING_MAKES_EMBEDDING_FUNCTIONAL",
        "MATERIAL_LOCKED_INVARIANCE_FAILED",
        "FIXED_WORLD_FIELD_EFFECT_WEAK",
    }


if __name__ == "__main__":
    main()
