import unittest

import numpy as np

from operaattori.cable_path import PassiveCableParams, path_abcd
from operaattori.order_nulls import (
    compose_orders,
    endpoint_preserving_permutations,
    full_permutations,
    monotone_radius_orders,
    precompute_segment_matrices,
    segment_midpoints_um,
    transfer_features,
    within_window_permutations,
)


class OrderNullTests(unittest.TestCase):
    def setUp(self):
        self.lengths = np.asarray([8., 12., 9., 15., 11., 7., 13., 10.])
        self.radii = np.asarray([1.8, 1.4, 1.5, 1.0, 0.9, 0.8, 0.65, 0.55])
        self.freqs = np.asarray([5., 40., 100.])

    def test_batch_composition_matches_reference(self):
        mats = precompute_segment_matrices(self.lengths, self.radii, self.freqs)
        orders = np.stack([
            np.arange(len(self.lengths)),
            np.arange(len(self.lengths))[::-1],
        ])
        got = compose_orders(mats, orders)
        for oi, order in enumerate(orders):
            for fi, f in enumerate(self.freqs):
                ref = path_abcd(self.lengths, self.radii, float(f), order=order)
                self.assertTrue(np.allclose(got[oi, fi], ref, rtol=1e-12, atol=1e-12))

    def test_full_permutations_preserve_multiset(self):
        rng = np.random.default_rng(1)
        p = full_permutations(12, 20, rng)
        target = np.arange(12)
        for row in p:
            self.assertTrue(np.array_equal(np.sort(row), target))

    def test_endpoint_shuffle_keeps_endpoint_zones(self):
        rng = np.random.default_rng(2)
        p = endpoint_preserving_permutations(self.lengths, 20, rng, fraction=0.20)
        mids = segment_midpoints_um(self.lengths)
        total = np.sum(self.lengths)
        fixed = (mids <= .2 * total) | (mids >= .8 * total)
        base = np.arange(len(self.lengths))
        self.assertTrue(np.all(p[:, fixed] == base[fixed]))

    def test_within_window_shuffle_never_crosses_window(self):
        rng = np.random.default_rng(3)
        p = within_window_permutations(self.lengths, 20.0, 20, rng)
        bins = np.floor(segment_midpoints_um(self.lengths) / 20.0).astype(int)
        for row in p:
            self.assertTrue(np.array_equal(bins[row], bins))

    def test_monotone_attackers_are_monotone(self):
        down, up = monotone_radius_orders(self.radii)
        self.assertTrue(np.all(np.diff(self.radii[down]) <= 0))
        self.assertTrue(np.all(np.diff(self.radii[up]) >= 0))

    def test_transfer_features_are_finite(self):
        mats = precompute_segment_matrices(self.lengths, self.radii, self.freqs)
        M = compose_orders(mats, np.arange(len(self.lengths)))
        feat = transfer_features(M)
        self.assertTrue(np.all(np.isfinite(feat)))


if __name__ == "__main__":
    unittest.main()
