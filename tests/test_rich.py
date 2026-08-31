import unittest

import numpy as np

from operaattori.rich import (
    BudgetMovingOperator,
    RichConfig,
    rich_drifting_world,
    run_seed,
)


class RichMovingOperatorTests(unittest.TestCase):
    def test_teacher_is_reproducible_and_finite(self):
        c = RichConfig(steps=400, burnin=40, seeds=1)
        a = rich_drifting_world(7, c)
        b = rich_drifting_world(7, c)
        for xa, xb in zip(a, b):
            self.assertTrue(np.allclose(xa, xb))
            self.assertTrue(np.all(np.isfinite(xa)))

    def test_moving_operator_stays_stable(self):
        c = RichConfig(steps=200, burnin=20, seeds=1)
        model = BudgetMovingOperator(4, c, move_operator=True)
        for t in range(200):
            p = model.predict(float(np.sin(t / 5)))
            model.learn(float(np.cos(t / 11)), p)
            lam = model.lambdas()
            self.assertTrue(np.all(lam > 0.0))
            self.assertTrue(np.all(lam < 1.0))

    def test_motion_helps_on_small_rich_world(self):
        c = RichConfig(
            steps=2600,
            burnin=300,
            seeds=1,
            moving_sizes=(2, 4),
            context_windows=(4, 8, 16),
        )
        row = run_seed(0, c)
        ratios = []
        for n in c.moving_sizes:
            m = row["moving"][str(n)]["mse"]
            f = row["frozen"][str(n)]["mse"]
            ratios.append(m / (f + 1e-12))
        self.assertLess(min(ratios), 1.0)

    def test_online_state_accounting(self):
        c = RichConfig(steps=100, burnin=10, seeds=1)
        self.assertEqual(
            BudgetMovingOperator(4, c, move_operator=True).online_scalars,
            20,
        )


if __name__ == "__main__":
    unittest.main()
