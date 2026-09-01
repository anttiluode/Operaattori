import unittest

import numpy as np

from operaattori.compartment_match import (
    greedy_dispersed_match,
    normalized_offdiagonal_coupling,
    select_compact_midpoint_sites,
    select_even_sites,
)


class CompartmentMatchTests(unittest.TestCase):
    def test_even_sites_are_unique_and_cover_ends(self):
        run = list(range(10, 110))
        s = select_even_sites(run, 12)
        self.assertEqual(len(s), 12)
        self.assertEqual(s[0], 10)
        self.assertEqual(s[-1], 109)
        self.assertEqual(len(np.unique(s)), len(s))

    def test_compact_midpoint_sites_respect_span(self):
        run = list(range(1, 21))
        lengths = np.ones(21) * 5.0
        sites, span = select_compact_midpoint_sites(run, lengths, 5, 30.0)
        self.assertEqual(len(sites), 5)
        self.assertLessEqual(span, 30.0 + 1e-12)
        self.assertGreater(sites[0], run[0])
        self.assertLess(sites[-1], run[-1])

    def test_match_avoids_target_run_and_spreads(self):
        n = 80
        run_ids = np.repeat(np.arange(8), 10)
        z = np.exp(np.linspace(3.0, 5.0, n))
        t = np.exp(np.linspace(-5.0, -1.0, n))
        target = np.asarray([0, 2, 4, 6, 8])
        pool = np.arange(n)
        chosen, diag = greedy_dispersed_match(
            target, 0, pool, run_ids, z, t, min_distinct_runs=4
        )
        self.assertTrue(np.all(run_ids[chosen] != 0))
        self.assertEqual(len(np.unique(chosen)), len(chosen))
        self.assertGreaterEqual(diag["distinct_match_runs"], 4)

    def test_offdiagonal_coupling_distinguishes_isolated(self):
        clustered = np.asarray([[100.0, 60.0], [60.0, 100.0]])
        dispersed = np.asarray([[100.0, 5.0], [5.0, 100.0]])
        self.assertGreater(
            normalized_offdiagonal_coupling(clustered),
            normalized_offdiagonal_coupling(dispersed),
        )


if __name__ == "__main__":
    unittest.main()
