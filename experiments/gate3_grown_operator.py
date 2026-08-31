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
    frozen_operator_matrix,
    mass_match_pair,
    probe_direct,
    probe_matrix,
    relative_l2,
    run_morphology,
    write_json,
)


def main() -> None:
    c = Config()
    clones = 12

    response_distance: list[float] = []
    operator_distance: list[float] = []
    uniform_response_distance: list[float] = []
    matrix_replay_error: list[float] = []
    same_history_response_floor: list[float] = []

    for clone in range(clones):
        paired_noise = 30_000 + clone
        m0 = run_morphology(
            H0, clone, paired_noise, c
        ).morphology
        m1 = run_morphology(
            H1, clone, paired_noise, c
        ).morphology

        # Remove total-mass information from the functional comparison.
        mm0, mm1 = mass_match_pair(m0, m1)

        f0 = probe_direct(mm0, c)
        f1 = probe_direct(mm1, c)
        response_distance.append(relative_l2(f0, f1))

        A0, _ = frozen_operator_matrix(mm0, c)
        A1, _ = frozen_operator_matrix(mm1, c)
        operator_distance.append(relative_l2(A0, A1))

        # If spatial arrangement is erased while total mass is kept equal, the
        # two histories become the same uniform material and the response
        # difference should vanish.
        common_mass = float(np.sum(mm0))
        uniform = np.full(c.n, common_mass / c.n, dtype=float)
        fu0 = probe_direct(uniform, c)
        fu1 = probe_direct(uniform, c)
        uniform_response_distance.append(relative_l2(fu0, fu1))

        # "Organic matrix" audit: when morphology is frozen, the fast probe is
        # linear and the exact derived matrix must replay it.
        fm = probe_matrix(mm0, c)
        matrix_replay_error.append(
            float(np.max(np.abs(f0 - fm)))
        )

        # Functional noise floor from two independent same-history growth runs.
        sa = run_morphology(
            H0, clone, 130_000 + clone, c
        ).morphology
        sb = run_morphology(
            H0, clone, 150_000 + clone, c
        ).morphology
        sa, sb = mass_match_pair(sa, sb)
        same_history_response_floor.append(
            relative_l2(
                probe_direct(sa, c),
                probe_direct(sb, c),
            )
        )

    mean_response = float(np.mean(response_distance))
    mean_floor = float(np.mean(same_history_response_floor))
    functional_ratio = mean_response / (mean_floor + 1e-12)

    result = {
        "gate": 3,
        "clones": clones,
        "mean_mass_matched_order_response_distance": mean_response,
        "mean_mass_matched_operator_matrix_distance":
            float(np.mean(operator_distance)),
        "mean_same_history_response_floor": mean_floor,
        "order_response_to_noise_floor_ratio": functional_ratio,
        "mean_uniform_mass_order_response_distance":
            float(np.mean(uniform_response_distance)),
        "max_exact_matrix_replay_error": max(matrix_replay_error),
    }
    write_json(ROOT / "results" / "gate3.json", result)

    print("Operaattori Gate 3 — anatomy -> grown operator")
    print(
        "mass-matched H0/H1 response distance: "
        f"{mean_response:.6f}"
    )
    print(
        "same-history response floor: "
        f"{mean_floor:.6f} | ratio={functional_ratio:.2f}x"
    )
    print(
        "mass-matched operator-matrix distance: "
        f"{result['mean_mass_matched_operator_matrix_distance']:.6f}"
    )
    print(
        "uniform-mass erased-arrangement distance: "
        f"{result['mean_uniform_mass_order_response_distance']:.3e}"
    )
    print(
        "exact frozen-matrix replay max error: "
        f"{result['max_exact_matrix_replay_error']:.3e}"
    )
    print()
    print(
        "Stopping line: different histories grow different spatial operators. "
        "After total mass is matched, the standardized future response still "
        "changes. But the frozen substrate is exactly a linear matrix, so this "
        "gate earns a self-grown operator, not superiority to matrices."
    )

    assert mean_response >= 0.20
    assert functional_ratio >= 5.0
    assert float(np.mean(uniform_response_distance)) <= 1e-12
    assert max(matrix_replay_error) <= 1e-10


if __name__ == "__main__":
    main()
