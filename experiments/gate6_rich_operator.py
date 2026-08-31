from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.rich import RichConfig, run_battery, write_json


def main() -> None:
    c = RichConfig()
    receipt = run_battery(c)
    write_json(ROOT / "results" / "gate6.json", receipt)

    s = receipt["summary"]

    print("Operaattori Gate 6 — out-of-family moving operator")
    print()
    print("Hidden world: six independently drifting memory modes")
    print()
    print(f"{'model':18s} {'online':>8s} {'MSE':>18s} {'horizon corr':>14s}")

    for n in c.moving_sizes:
        row = s["moving"][str(n)]
        print(
            f"{('moving-' + str(n)):18s} "
            f"{row['online_scalars']:8d} "
            f"{row['mse']['mean']:9.6g} ± {row['mse']['std']:.2g} "
            f"{row['horizon_correlation']['mean']:9.3f}"
        )
        f = s["frozen"][str(n)]
        print(
            f"{('frozen-' + str(n)):18s} "
            f"{f['online_scalars']:8d} "
            f"{f['mse']['mean']:9.6g} ± {f['mse']['std']:.2g} "
            f"{'-':>14s}"
        )

    for w in c.context_windows:
        row = s["context"][str(w)]
        print(
            f"{('context-' + str(w)):18s} "
            f"{row['online_scalars']:8d} "
            f"{row['mse']['mean']:9.6g} ± {row['mse']['std']:.2g} "
            f"{'-':>14s}"
        )

    print()
    print("Budget/error Pareto frontier:")
    for row in receipt["pareto_frontier"]:
        print(
            f"  {row['name']:16s} "
            f"online={row['online_scalars']:3d} "
            f"MSE={row['mse']:.6g}"
        )

    moving_improvements = []
    best_corr = -1.0
    for n in c.moving_sizes:
        m = s["moving"][str(n)]["mse"]["mean"]
        f = s["frozen"][str(n)]["mse"]["mean"]
        moving_improvements.append(m / (f + 1e-12))
        best_corr = max(
            best_corr,
            s["moving"][str(n)]["horizon_correlation"]["mean"],
        )

    frontier_names = {row["name"] for row in receipt["pareto_frontier"]}
    moving_on_frontier = any(
        f"moving-{n}" in frontier_names for n in c.moving_sizes
    )

    print()
    print(
        "Stopping line: Gate 6 asks whether moving-operator adaptation survives "
        "a teacher that is richer than the student's own family. A pass means "
        "motion still buys something over a frozen operator and at least one "
        "moving model lies on the online-state budget/error frontier. It does "
        "not mean the moving representation wins at every budget."
    )

    assert min(moving_improvements) <= 0.90
    assert best_corr >= 0.75
    assert moving_on_frontier


if __name__ == "__main__":
    main()
