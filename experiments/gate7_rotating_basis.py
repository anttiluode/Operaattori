from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.rotating import RotationConfig, run_battery, write_json


def main() -> None:
    c = RotationConfig()
    receipt = run_battery(c)
    write_json(ROOT / "results" / "gate7.json", receipt)

    s = receipt["summary"]

    print("Operaattori Gate 7 — move the basis")
    print()
    print("Hidden world: rotating two-mode operator")
    print(
        f"{'model':20s} {'online':>8s} {'MSE':>18s} "
        f"{'operator err':>14s} {'basis align':>12s}"
    )

    full = s["full"]
    diag = s["diagonal"]
    frozen = s["frozen"]

    print(
        f"{'moving-full':20s} {full['online_scalars']:8d} "
        f"{full['mse']['mean']:9.6g} ± {full['mse']['std']:.2g} "
        f"{full['operator_error']['mean']:9.4f} "
        f"{full['basis_alignment']['mean']:9.3f}"
    )
    print(
        f"{'moving-diagonal':20s} {diag['online_scalars']:8d} "
        f"{diag['mse']['mean']:9.6g} ± {diag['mse']['std']:.2g} "
        f"{diag['operator_error']['mean']:9.4f} "
        f"{'-':>12s}"
    )
    print(
        f"{'frozen-full':20s} {frozen['online_scalars']:8d} "
        f"{frozen['mse']['mean']:9.6g} ± {frozen['mse']['std']:.2g} "
        f"{'-':>14s} {'-':>12s}"
    )

    for w in c.context_windows:
        row = s["context"][str(w)]
        print(
            f"{('context-' + str(w)):20s} "
            f"{row['online_scalars']:8d} "
            f"{row['mse']['mean']:9.6g} ± {row['mse']['std']:.2g} "
            f"{'-':>14s} {'-':>12s}"
        )

    rls = s["rls32"]
    print(
        f"{'rls-32':20s} {rls['online_scalars']:8d} "
        f"{rls['mse']['mean']:9.6g} ± {rls['mse']['std']:.2g} "
        f"{'-':>14s} {'-':>12s}"
    )

    print()
    print("Budget/error Pareto frontier:")
    for row in receipt["pareto_frontier"]:
        print(
            f"  {row['name']:18s} "
            f"online={row['online_scalars']:4d} "
            f"MSE={row['mse']:.6g}"
        )

    full_mse = full["mse"]["mean"]
    diag_mse = diag["mse"]["mean"]
    frozen_mse = frozen["mse"]["mean"]

    print()
    print(
        "Stopping line: previous gates moved only eigenvalues. Gate 7 counts "
        "only if allowing the operator basis itself to move improves both "
        "prediction and direct matrix reconstruction relative to a moving "
        "diagonal attacker. Explicit-context models remain present as the "
        "non-operator alternative."
    )

    assert full_mse <= 0.80 * diag_mse
    assert full_mse <= 0.80 * frozen_mse
    assert (
        full["operator_error"]["mean"]
        <= 0.80 * diag["operator_error"]["mean"]
    )
    assert full["basis_alignment"]["mean"] >= 0.60


if __name__ == "__main__":
    main()
