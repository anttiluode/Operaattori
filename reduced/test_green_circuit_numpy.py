import unittest

import numpy as np

from reduced.green_circuit_numpy import (
    GreenCircuit,
    magnesium_block,
    shift_template,
    timed_conductance_matrix,
)


class ReducedGreenCircuitTests(unittest.TestCase):
    def test_shift_has_no_wraparound(self):
        x = np.asarray([1.0, 2.0, 3.0, 4.0])
        y = shift_template(x, 2)
        np.testing.assert_allclose(y, [0.0, 0.0, 1.0, 2.0])

    def test_timed_matrix(self):
        x = np.asarray([1.0, 2.0, 3.0])
        got = timed_conductance_matrix(x, [0, 1])
        np.testing.assert_allclose(
            got,
            [[1.0, 2.0, 3.0], [0.0, 1.0, 2.0]],
        )

    def test_magnesium_block_increases_with_depolarization(self):
        b = magnesium_block(np.asarray([-80.0, -50.0, -20.0]))
        self.assertTrue(np.all(np.diff(b) > 0.0))

    def test_zero_local_feedback_matches_direct_transport(self):
        n = 64
        local_h = np.zeros((1, 1, n))
        soma_h = np.zeros((1, n))
        soma_h[0, 0] = 2.0
        baseline = np.full((1, n), -70.0)

        circuit = GreenCircuit(
            local_h,
            soma_h,
            baseline,
            relative_tolerance=1e-12,
        )
        ga = np.zeros((1, n))
        ga[0, 4:8] = 0.001
        gn = np.zeros_like(ga)

        result = circuit.solve(ga, gn)
        expected_current = ga * 70.0
        expected_soma = 2.0 * expected_current[0]

        self.assertTrue(result.converged)
        np.testing.assert_allclose(
            result.current_nA,
            expected_current,
            rtol=1e-12,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            result.soma_depolarization_mV,
            expected_soma,
            rtol=1e-12,
            atol=1e-12,
        )

    def test_small_nonlinear_feedback_converges(self):
        n = 128
        local_h = np.zeros((2, 2, n))
        local_h[0, 0, 0] = 0.8
        local_h[1, 1, 0] = 0.8
        local_h[0, 1, 1] = 0.15
        local_h[1, 0, 1] = 0.15
        soma_h = np.zeros((2, n))
        soma_h[:, 0] = [0.5, 0.4]
        baseline = np.full((2, n), -70.0)

        circuit = GreenCircuit(local_h, soma_h, baseline)
        ga = np.zeros((2, n))
        gn = np.zeros((2, n))
        ga[:, 10:20] = 0.001
        gn[:, 10:50] = 0.002

        result = circuit.solve(ga, gn)
        self.assertTrue(result.converged)
        self.assertTrue(np.all(np.isfinite(result.current_nA)))
        self.assertGreater(
            float(np.max(result.soma_depolarization_mV)),
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
