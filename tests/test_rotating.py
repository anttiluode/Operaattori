import unittest

import numpy as np

from operaattori.rotating import (
    FullMovingOperator,
    RotationConfig,
    operator_from_coordinates,
    rotating_teacher,
    run_seed,
)


class RotatingOperatorTests(unittest.TestCase):
    def test_operator_stays_symmetric_stable(self):
        c = RotationConfig(steps=200, burnin=20, seeds=1)
        m = FullMovingOperator(c, move_basis=True, move_operator=True)
        for t in range(200):
            x = np.asarray([np.sin(t / 7), np.cos(t / 11)])
            target = np.asarray([np.cos(t / 9), np.sin(t / 13)])
            p = m.predict(x)
            m.learn(target, p)
            A = m.matrix()
            self.assertTrue(np.allclose(A, A.T, atol=1e-12))
            eig = np.linalg.eigvalsh(A)
            self.assertTrue(np.all(eig > 0.0))
            self.assertTrue(np.all(eig < 1.0))

    def test_angle_derivative_matches_finite_difference(self):
        logits = np.asarray([0.7, 2.0])
        phi = 0.43
        A, deriv = operator_from_coordinates(logits, phi)
        eps = 1e-6
        Ap, _ = operator_from_coordinates(logits, phi + eps)
        Am, _ = operator_from_coordinates(logits, phi - eps)
        numeric = (Ap - Am) / (2 * eps)
        self.assertTrue(np.allclose(deriv[2], numeric, atol=1e-6))
        self.assertEqual(A.shape, (2, 2))

    def test_teacher_reproducible(self):
        c = RotationConfig(steps=300, burnin=30, seeds=1)
        a = rotating_teacher(4, c)
        b = rotating_teacher(4, c)
        for key in ("x", "y", "phi", "A"):
            self.assertTrue(np.allclose(a[key], b[key]))

    def test_full_basis_can_beat_diagonal_on_small_world(self):
        c = RotationConfig(
            steps=2600,
            burnin=300,
            seeds=1,
            context_windows=(4, 8),
        )
        row = run_seed(0, c)
        self.assertLess(row["full"]["mse"], row["diagonal"]["mse"])


if __name__ == "__main__":
    unittest.main()
