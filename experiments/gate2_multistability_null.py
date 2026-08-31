from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.core import (
    Config,
    H0,
    H1,
    nearest_centroid_loo_pairs,
    run_morphology,
    unit_mass_shape,
    write_json,
)


def main() -> None:
    c = Config()
    clones = 24

    order_distances: list[float] = []
    same_history_distances: list[float] = []
    null0: list[np.ndarray] = []
    null1: list[np.ndarray] = []

    for clone in range(clones):
        # Paired order comparison: identical starting clone and the exact same
        # noise stream. Only history order differs.
        paired_noise = 20_000 + clone
        a = run_morphology(H0, clone, paired_noise, c).morphology
        b = run_morphology(H1, clone, paired_noise, c).morphology
        order_distances.append(
            float(
                np.linalg.norm(
                    unit_mass_shape(a) - unit_mass_shape(b)
                )
            )
        )

        # Multistability floor: same history, same initial clone, independently
        # resampled microscopic noise.
        a0 = run_morphology(H0, clone, 50_000 + clone, c).morphology
        a1 = run_morphology(H0, clone, 70_000 + clone, c).morphology
        b0 = run_morphology(H1, clone, 90_000 + clone, c).morphology
        b1 = run_morphology(H1, clone, 110_000 + clone, c).morphology

        same_history_distances.extend(
            [
                float(
                    np.linalg.norm(
                        unit_mass_shape(a0) - unit_mass_shape(a1)
                    )
                ),
                float(
                    np.linalg.norm(
                        unit_mass_shape(b0) - unit_mass_shape(b1)
                    )
                ),
            ]
        )

        # Same-history pseudo-classes: a classifier must not discover a stable
        # "class" from which independent noise stream was used.
        null0.append(unit_mass_shape(a0))
        null1.append(unit_mass_shape(a1))

    mean_order = float(np.mean(order_distances))
    mean_floor = float(np.mean(same_history_distances))
    ratio = mean_order / (mean_floor + 1e-12)
    null_acc = nearest_centroid_loo_pairs(
        np.asarray(null0), np.asarray(null1)
    )

    result = {
        "gate": 2,
        "clones": clones,
        "mean_paired_order_shape_distance": mean_order,
        "mean_same_history_multistability_distance": mean_floor,
        "order_to_multistability_ratio": ratio,
        "same_history_pseudoclass_acc": null_acc,
        "median_order_distance": float(np.median(order_distances)),
        "median_same_history_distance":
            float(np.median(same_history_distances)),
    }
    write_json(ROOT / "results" / "gate2.json", result)

    print("Operaattori Gate 2 — multistability noise floor")
    print(f"paired order shape distance: {mean_order:.6f}")
    print(f"same-history latch/noise floor: {mean_floor:.6f}")
    print(f"order / multistability floor: {ratio:.2f}x")
    print(f"same-history pseudo-class acc: {null_acc:.3f}")
    print()
    print(
        "Stopping line: the latch's own path dependence is not enough. "
        "The signal-order effect must clear the spread produced when the "
        "same history is rerun through independently noisy copies."
    )

    assert ratio >= 5.0
    assert 0.25 <= null_acc <= 0.75


if __name__ == "__main__":
    main()
