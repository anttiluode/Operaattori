from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.meta import MetaConfig, run_meta_battery, write_json


def main() -> None:
    c = MetaConfig()
    receipt = run_meta_battery(c)
    write_json(ROOT / "results" / "gate8.json", receipt)

    print("Operaattori Gate 8 — genes are the rule that moves A")
    print()
    print(
        f"outer-loop candidates: {receipt['candidate_count']} | "
        f"train worlds: {c.train_worlds} | test worlds: {c.test_worlds}"
    )
    print(f"selected candidate index: {receipt['selected_candidate_index']}")
    print(f"selected theta: {receipt['selected_theta']}")
    print()
    print(f"selected train MSE:          {receipt['selected_train_mse']:.6g}")
    print(
        "selected held-out test MSE: "
        f"{receipt['selected_test_mse']:.6g} ± "
        f"{receipt['selected_test_std']:.2g}"
    )
    print(f"hand-rule test MSE:          {receipt['hand_test_mse']:.6g}")
    print(f"frozen-operator test MSE:    {receipt['frozen_test_mse']:.6g}")
    print(
        "median random-rule test MSE:"
        f" {receipt['median_candidate_test_mse']:.6g}"
    )
    print(
        "per-world cheating oracle:   "
        f"{receipt['per_world_oracle_test_mse']:.6g}"
    )
    print(
        "64-token context LMS:        "
        f"{receipt['context64_test_mse']:.6g}"
    )
    print()
    print(
        "selected / oracle:           "
        f"{receipt['selected_to_oracle_ratio']:.3f}x"
    )
    print(
        "selected / hand:             "
        f"{receipt['selected_to_hand_ratio']:.3f}x"
    )
    print(
        "selected / median candidate: "
        f"{receipt['selected_to_median_candidate_ratio']:.3f}x"
    )
    print(
        "selected rule test rank:     "
        f"{receipt['selected_candidate_test_rank']} / "
        f"{receipt['candidate_count']}"
    )
    print(
        "held-out horizon correlation:"
        f" {receipt['selected_test_horizon_correlation']:.3f}"
    )
    print(
        "std of final operator tau across held-out worlds: "
        f"{receipt['test_final_operator_tau_std']:.3f}"
    )
    print()
    print(
        "Stopping line: theta is selected only from training worlds, then "
        "frozen. Every held-out world starts from the same A(0). A pass means "
        "one rule generates useful but different A(t) trajectories across "
        "unseen worlds. The per-world oracle is intentionally cheating and "
        "sets the upper bound on what this tiny rule family could do."
    )

    assert receipt["selected_to_median_candidate_ratio"] <= 0.80
    assert receipt["selected_to_hand_ratio"] <= 1.05
    assert receipt["selected_to_oracle_ratio"] <= 1.60
    assert receipt["selected_test_horizon_correlation"] >= 0.70
    assert receipt["selected_test_mse"] <= 0.50 * receipt["frozen_test_mse"]


if __name__ == "__main__":
    main()
