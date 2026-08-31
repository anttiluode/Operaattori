import unittest

import numpy as np

from operaattori.meta import (
    MetaConfig,
    meta_world,
    run_rule_on_world,
    sample_candidate_rules,
)


class MetaLawTests(unittest.TestCase):
    def test_world_is_reproducible(self):
        c = MetaConfig(
            steps=600,
            burnin=80,
            candidate_count=6,
            train_worlds=2,
            test_worlds=2,
        )
        a = meta_world(123, c)
        b = meta_world(123, c)
        for key in ("x", "y", "horizon"):
            self.assertTrue(np.allclose(a[key], b[key]))

    def test_candidate_population_is_reproducible(self):
        c = MetaConfig(candidate_count=8)
        self.assertEqual(
            sample_candidate_rules(c),
            sample_candidate_rules(c),
        )

    def test_same_rule_produces_different_operators_in_different_worlds(self):
        c = MetaConfig(
            steps=1200,
            burnin=150,
            candidate_count=6,
            train_worlds=2,
            test_worlds=2,
        )
        theta = sample_candidate_rules(c)[0]
        a = run_rule_on_world(theta, meta_world(1000, c), c)
        b = run_rule_on_world(theta, meta_world(1001, c), c)
        self.assertNotAlmostEqual(
            float(a["final_effective_tau"]),
            float(b["final_effective_tau"]),
            places=4,
        )

    def test_motion_beats_freeze_on_small_world(self):
        c = MetaConfig(
            steps=1600,
            burnin=200,
            candidate_count=6,
            train_worlds=2,
            test_worlds=2,
        )
        theta = sample_candidate_rules(c)[0]
        world = meta_world(1000, c)
        moving = run_rule_on_world(theta, world, c)
        frozen = run_rule_on_world(
            theta, world, c, move_operator=False
        )
        self.assertLess(moving["mse"], frozen["mse"])


if __name__ == "__main__":
    unittest.main()
