import unittest

import numpy as np

from operaattori.lie_path import (
    LiePathConfig,
    chain3_generators,
    commutator,
    lie_closure,
    ordered_flow,
    positive_loop,
    shear_generators,
)


class LiePathTests(unittest.TestCase):
    def test_closed_noncommuting_loop_leaves_residue(self):
        H, V = shear_generators()
        U = ordered_flow(positive_loop(H, V), 0.08)
        self.assertGreater(np.linalg.norm(U - np.eye(2)), 1e-3)
        self.assertLess(np.linalg.norm(H + V - H - V), 1e-12)

    def test_shear_commutator_is_saddle(self):
        H, V = shear_generators()
        expected = np.diag([1.0, -1.0])
        self.assertTrue(np.allclose(commutator(H, V), expected))

    def test_two_state_lie_closure_saturates_at_three(self):
        H, V = shear_generators()
        closure, _ = lie_closure([H, V])
        self.assertEqual(len(closure), 3)

    def test_nearest_neighbor_three_state_generators_fill_sl3(self):
        closure, _ = lie_closure(chain3_generators())
        self.assertEqual(len(closure), 8)

    def test_config_has_disjoint_train_and_test_loop_lengths(self):
        c = LiePathConfig()
        self.assertTrue(set(c.train_loop_counts).isdisjoint(c.test_loop_counts))


if __name__ == "__main__":
    unittest.main()
