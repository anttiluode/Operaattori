import unittest
import numpy as np

from operaattori.nmda_branch import (
    HUMAN, HUMAN_FROZEN_BLOCK, HUMAN_LINEAR_CURRENT,
    nmda_block, solve_equilibrium,
)


class NMDABranchTests(unittest.TestCase):
    def test_block_increases_with_depolarization(self):
        self.assertLess(float(nmda_block(-70.0, 0.078)), float(nmda_block(-30.0, 0.078)))

    def test_linear_current_obeys_superposition(self):
        Z = np.asarray([[100.0, 20.0], [20.0, 120.0]])
        one_a = solve_equilibrium(Z[:1,:1], np.asarray([1.0]), HUMAN_LINEAR_CURRENT)
        one_b = solve_equilibrium(Z[1:,1:], np.asarray([1.0]), HUMAN_LINEAR_CURRENT)
        both = solve_equilibrium(Z, np.ones(2), HUMAN_LINEAR_CURRENT)
        self.assertTrue(one_a['converged'])
        self.assertTrue(one_b['converged'])
        self.assertTrue(both['converged'])
        # Fixed currents are independent of voltage; only the passive Z mixes them.
        i_single = one_a['current_nA'][0] + one_b['current_nA'][0]
        self.assertAlmostEqual(float(np.sum(both['current_nA'])), float(i_single), places=12)

    def test_voltage_dependent_nmda_exceeds_frozen_block_at_high_resistance(self):
        Z = np.asarray([[250.0]])
        mult = np.asarray([20.0])
        full = solve_equilibrium(Z, mult, HUMAN)
        frozen = solve_equilibrium(Z, mult, HUMAN_FROZEN_BLOCK)
        self.assertTrue(full['converged'])
        self.assertTrue(frozen['converged'])
        self.assertGreater(full['depolarization_mV'][0], frozen['depolarization_mV'][0])


if __name__ == '__main__':
    unittest.main()
