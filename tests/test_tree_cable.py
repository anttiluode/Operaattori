import unittest

import numpy as np

from operaattori.tree_cable import (
    isolated_path_transfer,
    path_to_clamp,
    path_with_side_shunts_transfer,
    side_shunt_admittances,
    solve_tree_frequency,
)


class TreeCableTests(unittest.TestCase):
    def test_branchless_chain_matches_serial_abcd(self):
        # 0 is voltage-clamped root, 0->1->2 is a physical chain.
        parents = np.asarray([-1, 0, 1], dtype=np.int64)
        lengths = np.asarray([0.0, 80.0, 120.0])
        radii = np.asarray([1.0, 0.8, 0.6])
        active = np.asarray([True, True, True])
        clamped = np.asarray([True, False, False])

        state = solve_tree_frequency(
            parents, lengths, radii, active, clamped, 40.0
        )
        path = path_to_clamp(parents, clamped, 2)
        serial = isolated_path_transfer(path, lengths, radii, 40.0)
        self.assertAlmostEqual(abs(state.transfer_to_clamp[2] - serial), 0.0, places=11)

    def test_side_branch_changes_tip_transfer(self):
        # Main route 0->1->2 with side branch 1->3.
        parents = np.asarray([-1, 0, 1, 1], dtype=np.int64)
        lengths = np.asarray([0.0, 80.0, 120.0, 180.0])
        radii = np.asarray([1.0, 0.8, 0.6, 0.5])
        active = np.ones(4, dtype=bool)
        clamped = np.asarray([True, False, False, False])

        state = solve_tree_frequency(
            parents, lengths, radii, active, clamped, 80.0
        )
        path = path_to_clamp(parents, clamped, 2)
        isolated = isolated_path_transfer(path, lengths, radii, 80.0)
        self.assertGreater(
            abs(state.transfer_to_clamp[2] - isolated),
            1e-5,
        )

    def test_shunt_product_matches_whole_tree_elimination(self):
        parents = np.asarray([-1, 0, 1, 1, 3], dtype=np.int64)
        lengths = np.asarray([0.0, 70.0, 100.0, 90.0, 110.0])
        radii = np.asarray([1.0, 0.9, 0.6, 0.7, 0.5])
        active = np.ones(5, dtype=bool)
        clamped = np.asarray([True, False, False, False, False])
        f = 60.0

        state = solve_tree_frequency(parents, lengths, radii, active, clamped, f)
        path = path_to_clamp(parents, clamped, 4)
        shunts = side_shunt_admittances(
            parents, active, clamped, path, state
        )
        via_product = path_with_side_shunts_transfer(
            path, lengths, radii, shunts, f
        )
        self.assertAlmostEqual(
            abs(via_product - state.transfer_to_clamp[4]),
            0.0,
            places=10,
        )


if __name__ == "__main__":
    unittest.main()
