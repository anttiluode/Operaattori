from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.core import write_json
from operaattori.credit import CreditConfig, train


def mean_std(values):
    a = np.asarray(values, dtype=float)
    return {"mean": float(np.mean(a)), "std": float(np.std(a))}


def main() -> None:
    c = CreditConfig()
    seeds = tuple(range(12))
    arms = (
        "causal",
        "credit_shuffle",
        "eligibility_shuffle",
        "no_credit",
    )

    rows = {}
    for arm in arms:
        trials = [train(arm, seed, c) for seed in seeds]
        rows[arm] = {
            "selectivity": mean_std(
                [float(x["selectivity"]) for x in trials]
            ),
            "a_response": mean_std(
                [float(x["a_response"]) for x in trials]
            ),
            "b_response": mean_std(
                [float(x["b_response"]) for x in trials]
            ),
            "material_mass": mean_std(
                [float(np.sum(x["morphology"])) for x in trials]
            ),
        }

    causal = rows["causal"]["selectivity"]["mean"]
    shuffled_credit = rows["credit_shuffle"]["selectivity"]["mean"]
    shuffled_elig = rows["eligibility_shuffle"]["selectivity"]["mean"]
    no_credit = rows["no_credit"]["selectivity"]["mean"]

    result = {
        "gate": 4,
        "seeds": list(seeds),
        "config": c.__dict__,
        "arms": rows,
        "causal_minus_credit_shuffle": causal - shuffled_credit,
        "causal_minus_eligibility_shuffle": causal - shuffled_elig,
        "causal_minus_no_credit": causal - no_credit,
        "explicit_two_weight_digital_attacker_error": 0.0,
    }
    write_json(ROOT / "results" / "gate4.json", result)

    print("Operaattori Gate 4 — delayed consequence shapes operator")
    print(
        f"{'arm':24s} {'selectivity':>12s} "
        f"{'A peak':>12s} {'B peak':>12s} {'mass':>10s}"
    )
    for arm in arms:
        r = rows[arm]
        print(
            f"{arm:24s} "
            f"{r['selectivity']['mean']:8.3f}±"
            f"{r['selectivity']['std']:.3f} "
            f"{r['a_response']['mean']:12.6f} "
            f"{r['b_response']['mean']:12.6f} "
            f"{r['material_mass']['mean']:10.3f}"
        )

    print()
    print(
        "causal - shuffled consequence: "
        f"{causal - shuffled_credit:+.3f}"
    )
    print(
        "causal - shuffled eligibility: "
        f"{causal - shuffled_elig:+.3f}"
    )
    print(f"causal - no credit: {causal - no_credit:+.3f}")
    print()
    print(
        "Matrix ruler: an explicit two-weight digital mapping solves the "
        "A-high/B-low target exactly. Gate 4 is therefore only evidence that "
        "delayed scalar consequence can act through persistent local "
        "eligibility to grow a more selective physical operator."
    )

    assert causal >= 0.30
    assert causal - shuffled_credit >= 0.40
    assert causal - shuffled_elig >= 0.40
    assert causal - no_credit >= 0.60


if __name__ == "__main__":
    main()
