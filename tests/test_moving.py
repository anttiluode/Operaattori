import unittest

import numpy as np

from operaattori.moving import (
    MovingConfig,
    MovingDiagonalOperator,
    drifting_delay_world,
    run_delay_seed,
    run_kernel_seed,
)


class MovingOperatorTests(unittest.TestCase):
    def test_operator_is_stable_diagonal_matrix(self):
        c = MovingConfig(steps=64, burnin=8, seeds=1)
        m = MovingDiagonalOperator(c)
        for t in range(64):
            p = m.predict(float(np.sin(t)))
            m.learn(float(np.cos(t / 3)), p)
            eig = np.diag(m.operator_matrix())
            self.assertTrue(np.all(eig > 0.0))
            self.assertTrue(np.all(eig < 1.0))

    def test_small_kernel_world_rewards_operator_motion(self):
        c = MovingConfig(
            steps=2500,
            burnin=300,
            seeds=1,
            context_kernel=48,
            context_delay=24,
        )
        row = run_kernel_seed(0, c)
        self.assertLess(row["moving_mse"], row["frozen_mse"])
        self.assertGreater(row["operator_timescale_correlation"], 0.95)

    def test_exact_delay_is_a_real_attention_attack(self):
        c = MovingConfig(
            steps=2500,
            burnin=300,
            seeds=1,
            context_kernel=48,
            context_delay=24,
        )
        row = run_delay_seed(0, c)
        self.assertLess(row["attention_mse"], row["moving_mse"])
        self.assertGreater(row["attention_lag_correlation"], 0.95)

    def test_delay_teacher_never_exposes_future_input(self):
        c = MovingConfig(steps=200, burnin=20, seeds=1, context_delay=24)
        x, y, lag = drifting_delay_world(3, c)
        self.assertTrue(np.all(lag >= 1))
        self.assertEqual(len(x), len(y))
        self.assertEqual(len(y), c.steps)


if __name__ == "__main__":
    unittest.main()
