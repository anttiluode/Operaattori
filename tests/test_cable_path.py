import unittest

import numpy as np

from operaattori.cable_path import (
    PassiveCableParams,
    area_preserving_radius,
    cable_abcd,
    commutator_action_score,
    path_abcd,
    relative_complex_difference,
    sealed_distal_voltage_gain,
    sealed_input_impedance,
)


class CablePathTests(unittest.TestCase):
    def test_uniform_cable_segments_commute(self):
        M1 = cable_abcd(80.0, 0.8, 40.0)
        M2 = cable_abcd(150.0, 0.8, 40.0)
        self.assertLess(np.linalg.norm(M1 @ M2 - M2 @ M1), 1e-8)
        self.assertLess(commutator_action_score(M1, M2), 1e-12)

    def test_different_radii_are_noncommuting(self):
        M1 = cable_abcd(120.0, 0.35, 50.0)
        M2 = cable_abcd(120.0, 1.8, 50.0)
        self.assertGreater(commutator_action_score(M1, M2), 1e-5)

    def test_order_changes_function_with_same_segment_multiset(self):
        lengths = np.asarray([90.0, 140.0, 110.0, 170.0])
        radii = np.asarray([0.35, 1.2, 0.5, 1.8])
        f = 80.0
        forward = path_abcd(lengths, radii, f)
        reverse = path_abcd(lengths, radii, f, order=np.arange(len(lengths))[::-1])
        dz = relative_complex_difference(
            sealed_input_impedance(forward),
            sealed_input_impedance(reverse),
        )
        dg = relative_complex_difference(
            sealed_distal_voltage_gain(forward),
            sealed_distal_voltage_gain(reverse),
        )
        self.assertGreater(max(dz, dg), 1e-4)

    def test_area_matched_uniform_control_is_order_invariant(self):
        lengths = np.asarray([90.0, 140.0, 110.0, 170.0])
        radii = np.asarray([0.35, 1.2, 0.5, 1.8])
        r = area_preserving_radius(lengths, radii)
        rr = np.full_like(radii, r)
        f = 80.0
        forward = path_abcd(lengths, rr, f)
        reverse = path_abcd(lengths, rr, f, order=np.arange(len(lengths))[::-1])
        self.assertLess(
            relative_complex_difference(
                sealed_input_impedance(forward),
                sealed_input_impedance(reverse),
            ),
            1e-12,
        )


if __name__ == "__main__":
    unittest.main()
