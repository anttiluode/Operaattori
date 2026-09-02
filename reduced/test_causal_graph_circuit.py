from __future__ import annotations

import unittest

import numpy as np

from reduced.causal_graph_circuit import (
    CausalGraphCircuit,
    compile_dense_passive_graph,
    current_derivative_nA_per_mV,
    inward_synaptic_current_nA,
)


def tiny_graph():
    g_leak = 0.02
    g_ax = 0.12
    G = np.array(
        [
            [g_leak + g_ax, -g_ax, 0.0, 0.0],
            [-g_ax, g_leak + 2 * g_ax, -g_ax, 0.0],
            [0.0, -g_ax, g_leak + 2 * g_ax, -g_ax],
            [0.0, 0.0, -g_ax, g_leak + g_ax],
        ],
        dtype=float,
    )
    C = np.array([0.08, 0.07, 0.06, 0.08], dtype=float)
    return G, C


def pulse(t, onset, tau_r=0.3, tau_d=5.0):
    x = t - float(onset)
    y = np.where(
        x > 0.0,
        np.exp(-x / tau_d) - np.exp(-x / tau_r),
        0.0,
    )
    m = float(np.max(y))
    return y / m if m > 0 else y


def conductance_program(ntime=600, dt=0.05):
    t = np.arange(ntime, dtype=float) * dt
    ga = np.stack(
        [
            0.0010 * pulse(t, 4.0),
            0.0010 * pulse(t, 5.0),
            0.0010 * pulse(t, 6.0),
        ],
        axis=0,
    )
    gn = np.stack(
        [
            0.0030 * pulse(t, 4.0),
            0.0030 * pulse(t, 5.0),
            0.0030 * pulse(t, 6.0),
        ],
        axis=0,
    )
    return ga, gn


class CausalGraphCircuitTests(unittest.TestCase):
    def test_dense_compiler_matches_backward_euler_step(self):
        G, C = tiny_graph()
        sites = np.array([1, 2, 3], dtype=int)
        P, X = compile_dense_passive_graph(
            G, C, sites, dt_ms=0.05
        )

        state = np.array([0.1, 0.2, 0.3, 0.4])
        current = np.array([0.01, 0.02, 0.03])

        d = C / 0.05
        A = G + np.diag(d)
        B = np.zeros((4, 3), dtype=float)
        B[sites, np.arange(3)] = 1.0
        direct = np.linalg.solve(
            A, d * state + B @ current
        )
        compiled = P @ state + X @ current
        self.assertTrue(
            np.allclose(direct, compiled, rtol=0, atol=1e-14)
        )

    def test_synapse_derivative_matches_finite_difference(self):
        v = np.array([-72.0, -61.0, -47.0])
        ga = np.array([0.0010, 0.0014, 0.0008])
        gn = np.array([0.0030, 0.0025, 0.0035])
        analytic = current_derivative_nA_per_mV(
            v, ga, gn
        )
        eps = 1e-6
        numeric = (
            inward_synaptic_current_nA(v + eps, ga, gn)
            - inward_synaptic_current_nA(v - eps, ga, gn)
        ) / (2.0 * eps)
        self.assertTrue(
            np.allclose(analytic, numeric, rtol=2e-7, atol=2e-10)
        )

    def test_zero_conductance_stays_at_rest(self):
        G, C = tiny_graph()
        circuit = CausalGraphCircuit.from_dense_passive_graph(
            G, C, [1, 2, 3], 0, dt_ms=0.05
        )
        result = circuit.run(
            np.zeros((3, 20)),
            np.zeros((3, 20)),
        )
        self.assertTrue(result.all_steps_converged)
        self.assertEqual(result.max_newton_iterations, 0)
        self.assertEqual(
            float(np.max(np.abs(result.soma_depolarization_mV))),
            0.0,
        )

    def test_nonlinear_program_converges_and_closes_site_state(self):
        G, C = tiny_graph()
        circuit = CausalGraphCircuit.from_dense_passive_graph(
            G, C, [1, 2, 3], 0, dt_ms=0.05
        )
        ga, gn = conductance_program()
        result = circuit.run(ga, gn)

        self.assertTrue(result.all_steps_converged)
        self.assertLessEqual(result.max_newton_iterations, 4)
        self.assertEqual(result.total_line_search_failures, 0)
        self.assertLess(result.max_site_consistency_mV, 2e-10)
        self.assertGreater(
            float(np.max(result.soma_depolarization_mV)),
            0.5,
        )

    def test_compartment_relabeling_is_equivariant(self):
        G, C = tiny_graph()
        sites = np.array([1, 2, 3], dtype=int)
        soma = 0
        ga, gn = conductance_program(ntime=300)

        base = CausalGraphCircuit.from_dense_passive_graph(
            G, C, sites, soma, dt_ms=0.05
        ).run(ga, gn)

        # new-index -> old-index permutation
        perm = np.array([2, 0, 3, 1], dtype=int)
        old_to_new = np.argsort(perm)
        Gp = G[np.ix_(perm, perm)]
        Cp = C[perm]
        sites_p = old_to_new[sites]
        soma_p = int(old_to_new[soma])

        relabeled = CausalGraphCircuit.from_dense_passive_graph(
            Gp, Cp, sites_p, soma_p, dt_ms=0.05
        ).run(ga, gn)

        self.assertTrue(
            np.allclose(
                base.soma_depolarization_mV,
                relabeled.soma_depolarization_mV,
                rtol=1e-12,
                atol=1e-12,
            )
        )
        self.assertTrue(
            np.allclose(
                base.current_nA,
                relabeled.current_nA,
                rtol=1e-12,
                atol=1e-12,
            )
        )


if __name__ == "__main__":
    unittest.main()
