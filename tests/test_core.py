import unittest

import numpy as np

from operaattori.core import (
    Config,
    H0,
    H1,
    frozen_operator_matrix,
    history_invariants,
    mass_match_pair,
    probe_direct,
    probe_matrix,
    schedule,
)


class ProtocolTests(unittest.TestCase):
    def test_histories_match_marginals_first_and_last(self):
        inv = history_invariants()
        self.assertTrue(all(inv.values()))

    def test_complete_suffix_is_identical(self):
        c = Config()
        p0 = schedule(H0, c)
        p1 = schedule(H1, c)
        suffix = c.common_suffix_steps + c.washout_steps
        self.assertTrue(np.array_equal(p0[-suffix:], p1[-suffix:]))

    def test_washout_exceeds_seven_eligibility_time_constants(self):
        self.assertGreaterEqual(
            Config().washout_eligibility_time_constants, 7.0
        )


class OperatorTests(unittest.TestCase):
    def test_mass_match_really_matches_total_material(self):
        a = np.asarray([0.2, 0.4, 0.6])
        b = np.asarray([0.1, 0.1, 0.2])
        aa, bb = mass_match_pair(a, b)
        self.assertAlmostEqual(float(np.sum(aa)), float(np.sum(bb)))

    def test_frozen_operator_is_stable_enough_for_probe(self):
        c = Config()
        m = np.full(c.n, 0.5)
        A, _ = frozen_operator_matrix(m, c)
        rho = float(np.max(np.abs(np.linalg.eigvals(A))))
        self.assertLessEqual(rho, 1.0 + 1e-10)

    def test_matrix_replay_is_exact(self):
        c = Config()
        m = np.linspace(0.0, 0.9, c.n)
        direct = probe_direct(m, c)
        matrix = probe_matrix(m, c)
        self.assertLess(
            float(np.max(np.abs(direct - matrix))),
            1e-10,
        )


if __name__ == "__main__":
    unittest.main()
