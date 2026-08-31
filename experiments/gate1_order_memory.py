from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.core import (
    BistableStateRuler,
    Config,
    ContractiveMatrixRuler,
    H0,
    H1,
    nearest_centroid_loo_pairs,
    run_morphology,
    schedule,
    unit_mass_shape,
    write_json,
)


def ruler_accuracy(cls, config: Config, clones: int) -> float:
    arms = [[], []]
    for clone in range(clones):
        noise_seed = 100_000 + clone
        for arm, history in enumerate((H0, H1)):
            arms[arm].append(
                cls(config, clone, noise_seed).run(
                    schedule(history, config)
                )
            )
    return nearest_centroid_loo_pairs(
        np.asarray(arms[0]), np.asarray(arms[1])
    )


def main() -> None:
    c = Config()
    clones = 24
    raw = [[], []]
    shape = [[], []]
    fast_residual = []
    eligibility_residual = []

    for clone in range(clones):
        # Identical clone + identical microscopic noise schedule. The only
        # paired difference is the order protocol.
        noise_seed = 100_000 + clone
        for arm, history in enumerate((H0, H1)):
            world = run_morphology(history, clone, noise_seed, c)
            raw[arm].append(world.morphology.copy())
            shape[arm].append(unit_mass_shape(world.morphology))
            fast_residual.append(float(np.max(np.abs(world.fast))))
            eligibility_residual.append(
                float(np.max(np.abs(world.eligibility)))
            )

    raw0, raw1 = map(np.asarray, raw)
    shape0, shape1 = map(np.asarray, shape)
    raw_acc = nearest_centroid_loo_pairs(raw0, raw1)
    shape_acc = nearest_centroid_loo_pairs(shape0, shape1)

    contractive_acc = ruler_accuracy(
        ContractiveMatrixRuler, c, clones
    )
    bistable_acc = ruler_accuracy(BistableStateRuler, c, clones)

    result = {
        "gate": 1,
        "clones": clones,
        "raw_anatomy_nearest_centroid_acc": raw_acc,
        "unit_mass_shape_nearest_centroid_acc": shape_acc,
        "max_fast_residual": max(fast_residual),
        "max_eligibility_residual": max(eligibility_residual),
        "washout_eligibility_time_constants":
            c.washout_eligibility_time_constants,
        "contractive_matrix_ruler_acc": contractive_acc,
        "bistable_state_ruler_acc": bistable_acc,
        "mean_total_mass_h0": float(np.mean(np.sum(raw0, axis=1))),
        "mean_total_mass_h1": float(np.mean(np.sum(raw1, axis=1))),
    }
    write_json(ROOT / "results" / "gate1.json", result)

    print("Operaattori Gate 1 — signal -> persistent morphology")
    print(f"raw anatomy LOO nearest-centroid: {raw_acc:.3f}")
    print(f"unit-mass shape LOO nearest-centroid: {shape_acc:.3f}")
    print(
        "fast residual: "
        f"{max(fast_residual):.3e} | eligibility residual: "
        f"{max(eligibility_residual):.3e}"
    )
    print(
        "mean total mass H0/H1: "
        f"{result['mean_total_mass_h0']:.3f} / "
        f"{result['mean_total_mass_h1']:.3f}"
    )
    print(
        "rulers — contractive matrix: "
        f"{contractive_acc:.3f}, bistable abstract state: "
        f"{bistable_acc:.3f}"
    )
    print()
    print(
        "Stopping line: order is decodable from persistent spatial shape even "
        "after total mass is normalized. This is morphological memory. The "
        "bistable abstract ruler also remembers, so no unique computation "
        "claim is earned."
    )

    assert raw_acc >= 0.90
    assert shape_acc >= 0.90
    assert max(fast_residual) <= 1e-5
    assert max(eligibility_residual) <= 0.005
    assert contractive_acc <= 0.65


if __name__ == "__main__":
    main()
