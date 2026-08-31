from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.moving import MovingConfig, run_battery, write_json


def pm(summary, key):
    row = summary[key]
    return f"{row['mean']:.6g} ± {row['std']:.2g}"


def main() -> None:
    c = MovingConfig()
    receipt = run_battery(c)
    write_json(ROOT / "results" / "gate5.json", receipt)

    k = receipt["kernel"]["summary"]
    d = receipt["delay_attack"]["summary"]

    print("Operaattori Gate 5 — the matrix itself is state")
    print()
    print("WORLD A: smoothly drifting memory kernel")
    print(f"moving operator MSE:          {pm(k, 'moving_mse')}")
    print(f"frozen same-state SSM MSE:   {pm(k, 'frozen_mse')}")
    print(f"64-token context LMS MSE:    {pm(k, 'explicit_context_mse')}")
    print(
        "moving / frozen MSE:        "
        f"{k['moving_vs_frozen_ratio']['mean']:.3f}x"
    )
    print(
        "moving / context MSE:       "
        f"{k['moving_vs_context_ratio']['mean']:.3f}x"
    )
    print(
        "operator-timescale / world correlation: "
        f"{k['operator_timescale_correlation']['mean']:.3f}"
    )
    print(
        "online stored/mutable scalars: moving=10, "
        "explicit-context=128"
    )
    print()

    print("WORLD B ATTACK: exact drifting delay")
    print(f"moving operator MSE:          {pm(d, 'moving_mse')}")
    print(f"frozen same-state SSM MSE:   {pm(d, 'frozen_mse')}")
    print(f"32-token context LMS MSE:    {pm(d, 'explicit_context_mse')}")
    print(f"lag-attention MSE:           {pm(d, 'attention_mse')}")
    print(
        "moving effective-timescale / true-lag correlation: "
        f"{d['moving_lag_correlation']['mean']:.3f}"
    )
    print(
        "attention estimated-lag / true-lag correlation:     "
        f"{d['attention_lag_correlation']['mean']:.3f}"
    )
    print()
    print(
        "Stopping line: a tiny moving recurrent operator earns a real niche "
        "on a smoothly drifting operator family and compresses history into "
        "10 online scalars. It does NOT replace explicit context: on exact "
        "moving delays, attention wins decisively and the moving operator "
        "does not track the useful lag."
    )

    assert k["moving_vs_frozen_ratio"]["mean"] <= 0.50
    assert k["moving_vs_context_ratio"]["mean"] <= 0.97
    assert k["operator_timescale_correlation"]["mean"] >= 0.98

    assert d["attention_mse"]["mean"] <= 0.20 * d["moving_mse"]["mean"]
    assert d["explicit_context_mse"]["mean"] <= 0.35 * d["moving_mse"]["mean"]
    assert d["attention_lag_correlation"]["mean"] >= 0.98


if __name__ == "__main__":
    main()
