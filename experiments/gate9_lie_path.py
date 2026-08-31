from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.lie_path import LiePathConfig, run_battery, write_json


def main() -> None:
    c = LiePathConfig()
    receipt = run_battery(c)
    write_json(ROOT / "results" / "gate9.json", receipt)

    loop = receipt["closed_loop"]
    dec = receipt["held_out_path_decoder"]
    lie = receipt["lie_closure"]

    print("Operaattori Gate 9 — path through operator space")
    print()
    print("Closed noncommuting loop: H -> V -> -H -> -V")
    print(f"net integrated generator:            {loop['net_exposure_norm']:.3e}")
    print(f"noncommuting loop residue:           {loop['positive_residue_norm']:.6g}")
    print(f"commuting-control loop residue:      {loop['commuting_positive_residue_norm']:.3e}")
    print(
        "reverse-loop delta / [H,V] alignment: "
        f"{loop['reverse_difference_commutator_alignment']:.6f}"
    )
    print(
        "residue scaling exponent vs epsilon: "
        f"{loop['epsilon_scaling_slope']:.4f}"
    )
    print()
    print("Held-out path-length decoder")
    print(
        "noncommuting final-state corr:        "
        f"{dec['correlation']['mean']:.4f} ± {dec['correlation']['std']:.4f}"
    )
    print(
        "noncommuting final-state R^2:         "
        f"{dec['r2']['mean']:.4f} ± {dec['r2']['std']:.4f}"
    )
    print(
        "commuting-control final-state R^2:    "
        f"{dec['commuting_r2']['mean']:.4f} ± {dec['commuting_r2']['std']:.4f}"
    )
    print("digital counter attacker:             exact by construction")
    print()
    print("Lie closure")
    print(
        "2-state shear primitives -> closure: "
        f"{lie['two_state_primitive_dimension']} -> {lie['two_state_closure_dimension']}"
    )
    print(
        "3-state nearest-neighbor primitives -> closure: "
        f"{lie['three_state_primitive_dimension']} -> {lie['three_state_closure_dimension']} "
        f"(full sl(3) = {lie['full_sl3_dimension']})"
    )
    print(
        "non-neighbor E13 distance to primitive span: "
        f"{lie['E13_distance_to_primitive_span']:.3g}"
    )
    print(
        "non-neighbor E13 distance to closure span:   "
        f"{lie['E13_distance_to_closure_span']:.3g}"
    )
    print()
    print(
        "Stopping line: order can live in the time-ordered flow even after the "
        "instantaneous generator returns to neutral and its first-order integral "
        "is zero. The positive mechanism is noncommutativity / Lie brackets. "
        "The negative boundary is equally important: a fixed 2x2 linear system "
        "closes after three traceless directions, and a tiny digital counter "
        "solves the toy path target exactly. This is a factorization/development "
        "mechanism, not a new function class."
    )

    assert loop["net_exposure_norm"] <= 1e-12
    assert loop["positive_residue_norm"] >= 1e-3
    assert loop["commuting_positive_residue_norm"] <= 1e-10
    assert loop["reverse_difference_commutator_alignment"] >= 0.98
    assert 1.90 <= loop["epsilon_scaling_slope"] <= 2.10
    assert dec["correlation"]["mean"] >= 0.94
    assert dec["r2"]["mean"] >= 0.85
    assert dec["commuting_r2"]["mean"] <= 0.05
    assert lie["two_state_closure_dimension"] == 3
    assert lie["three_state_closure_dimension"] == 8
    assert lie["E13_distance_to_primitive_span"] >= 0.90
    assert lie["E13_distance_to_closure_span"] <= 1e-10


if __name__ == "__main__":
    main()
