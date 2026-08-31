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
    history_invariants,
    nearest_centroid_loo_pairs,
    schedule,
    write_json,
)


def collect(config: Config, clones: int = 24):
    contractive = [[], []]
    bistable = [[], []]
    for clone in range(clones):
        for arm, history in enumerate((H0, H1)):
            noise_seed = 100_000 + clone  # exact paired noise schedule
            protocol = schedule(history, config)
            contractive[arm].append(
                ContractiveMatrixRuler(
                    config, clone, noise_seed
                ).run(protocol)
            )
            bistable[arm].append(
                BistableStateRuler(
                    config, clone, noise_seed
                ).run(protocol)
            )
    return (
        np.asarray(contractive[0]),
        np.asarray(contractive[1]),
        np.asarray(bistable[0]),
        np.asarray(bistable[1]),
    )


def main() -> None:
    c = Config()
    inv = history_invariants()
    p0 = schedule(H0, c)
    p1 = schedule(H1, c)

    suffix = c.common_suffix_steps + c.washout_steps
    same_suffix = bool(np.array_equal(p0[-suffix:], p1[-suffix:]))

    x0, x1, b0, b1 = collect(c)
    contractive_acc = nearest_centroid_loo_pairs(x0, x1)
    bistable_acc = nearest_centroid_loo_pairs(b0, b1)

    result = {
        "gate": 0,
        "histories": [H0, H1],
        "invariants": inv,
        "same_complete_suffix": same_suffix,
        "washout_eligibility_time_constants":
            c.washout_eligibility_time_constants,
        "contractive_matrix_nearest_centroid_acc": contractive_acc,
        "bistable_state_nearest_centroid_acc": bistable_acc,
    }
    write_json(ROOT / "results" / "gate0.json", result)

    print("Operaattori Gate 0 — protocol + rulers")
    print(f"histories: {H0} vs {H1}")
    print(f"invariants: {inv}")
    print(f"exact common suffix: {same_suffix}")
    print(
        "washout / eligibility tau: "
        f"{c.washout_eligibility_time_constants:.2f}"
    )
    print(f"contractive matrix ruler acc: {contractive_acc:.3f}")
    print(f"bistable state ruler acc:     {bistable_acc:.3f}")
    print()
    print(
        "Interpretation: contraction should forget after the common suffix. "
        "A same-capacity hysteretic abstract state is allowed to remember. "
        "Morphology must therefore earn more than 'hysteresis stores order'."
    )

    assert all(bool(v) for v in inv.values())
    assert same_suffix
    assert c.washout_eligibility_time_constants >= 7.0
    assert contractive_acc <= 0.65
    # This dangerous null should be able to pass the storage assay.
    assert bistable_acc >= 0.90


if __name__ == "__main__":
    main()
