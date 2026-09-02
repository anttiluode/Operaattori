"""Pure-NumPy causal circuit runtime for a compiled passive cable graph.

Operaattori's passing cross-cell audit used a sparse backward-Euler cable solve
plus a local implicit NMDA solve.  This module extracts the simulator-free
runtime boundary:

    passive state_(n+1) = P @ state_n
    site voltage        = S passive + R @ J(site voltage)
    state_(n+1)         = passive + X @ J
    soma output         = state_(n+1)[soma_node]

where P and X are compiled once from the passive graph and
R = X[site_nodes, :].

The runtime itself depends only on NumPy.  ``from_dense_passive_graph`` is a
small reference compiler for modest graphs.  The 24-cell scientific audit used
sparse LU during compilation; large morphologies should keep that sparse
compiler rather than materializing a dense P matrix.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def magnesium_block(
    voltage_mV: np.ndarray,
    *,
    gamma_per_mV: float = 0.078,
    mg_divisor: float = 3.57,
) -> np.ndarray:
    """Released-model Jahr-Stevens-style magnesium block."""
    v = np.asarray(voltage_mV, dtype=float)
    return 1.0 / (
        1.0 + np.exp(-float(gamma_per_mV) * v) / float(mg_divisor)
    )


def inward_synaptic_current_nA(
    voltage_mV: np.ndarray,
    g_ampa_uS: np.ndarray,
    g_nmda_raw_uS: np.ndarray,
    *,
    reversal_mV: float = 0.0,
    gamma_per_mV: float = 0.078,
    mg_divisor: float = 3.57,
) -> np.ndarray:
    """Positive-inward AMPA+NMDA current in nA (uS*mV == nA)."""
    v = np.asarray(voltage_mV, dtype=float)
    ga = np.asarray(g_ampa_uS, dtype=float)
    gn = np.asarray(g_nmda_raw_uS, dtype=float)
    b = magnesium_block(
        v,
        gamma_per_mV=gamma_per_mV,
        mg_divisor=mg_divisor,
    )
    return (ga + gn * b) * (float(reversal_mV) - v)


def current_derivative_nA_per_mV(
    voltage_mV: np.ndarray,
    g_ampa_uS: np.ndarray,
    g_nmda_raw_uS: np.ndarray,
    *,
    reversal_mV: float = 0.0,
    gamma_per_mV: float = 0.078,
    mg_divisor: float = 3.57,
) -> np.ndarray:
    """Analytic derivative dJ/dV for the local Newton solve."""
    v = np.asarray(voltage_mV, dtype=float)
    ga = np.asarray(g_ampa_uS, dtype=float)
    gn = np.asarray(g_nmda_raw_uS, dtype=float)
    b = magnesium_block(
        v,
        gamma_per_mV=gamma_per_mV,
        mg_divisor=mg_divisor,
    )
    bp = float(gamma_per_mV) * b * (1.0 - b)
    return (
        gn * bp * (float(reversal_mV) - v)
        - (ga + gn * b)
    )


def compile_dense_passive_graph(
    G_uS: np.ndarray,
    C_nF: np.ndarray,
    site_nodes: list[int] | tuple[int, ...] | np.ndarray,
    *,
    dt_ms: float = 0.025,
) -> tuple[np.ndarray, np.ndarray]:
    """Compile a modest dense passive graph into one-step matrices P and X.

    The passive graph is

        C dv/dt + G v = I

    and backward Euler gives

        A = diag(C/dt) + G
        P = A^-1 diag(C/dt)
        X = A^-1 B

    with B injecting unit current at ``site_nodes``.

    Notes
    -----
    This helper intentionally uses only NumPy.  It materializes P as an n x n
    matrix and is therefore a reference/demo compiler, not the memory-efficient
    path used by Operaattori's large-morphology audit.
    """
    if float(dt_ms) <= 0:
        raise ValueError("dt_ms must be positive")
    G = np.asarray(G_uS, dtype=float)
    C = np.asarray(C_nF, dtype=float)
    if C.ndim != 1 or len(C) == 0:
        raise ValueError("C_nF must be a non-empty 1-D array")
    n = len(C)
    if G.shape != (n, n):
        raise ValueError(f"G_uS must have shape {(n, n)}")
    if not np.all(np.isfinite(G)) or not np.all(np.isfinite(C)):
        raise FloatingPointError("passive graph contains non-finite values")
    if np.any(C <= 0):
        raise ValueError("C_nF must be strictly positive")

    sites = np.asarray(site_nodes, dtype=int)
    if sites.ndim != 1 or len(sites) == 0:
        raise ValueError("site_nodes must be a non-empty 1-D sequence")
    if np.any(sites < 0) or np.any(sites >= n):
        raise IndexError("site_nodes contain an invalid compartment index")
    if len(np.unique(sites)) != len(sites):
        raise ValueError("site_nodes must be unique")

    d = C / float(dt_ms)
    A = G + np.diag(d)
    B = np.zeros((n, len(sites)), dtype=float)
    B[sites, np.arange(len(sites))] = 1.0

    rhs = np.concatenate((np.diag(d), B), axis=1)
    solved = np.linalg.solve(A, rhs)
    return solved[:, :n], solved[:, n:]


@dataclass(frozen=True)
class StepSolve:
    site_depolarization_mV: np.ndarray
    current_nA: np.ndarray
    converged: bool
    iterations: int
    residual_inf_mV: float
    line_search_failures: int
    failure: str | None = None


@dataclass(frozen=True)
class CausalSolveResult:
    current_nA: np.ndarray
    local_voltage_mV: np.ndarray
    soma_depolarization_mV: np.ndarray
    final_state_depolarization_mV: np.ndarray
    all_steps_converged: bool
    max_newton_iterations: int
    max_residual_inf_mV: float
    total_line_search_failures: int
    max_site_consistency_mV: float
    first_failure_time_index: int | None


class CausalGraphCircuit:
    """Run a compiled passive state-space operator with local NMDA feedback."""

    def __init__(
        self,
        passive_step_matrix: np.ndarray,
        input_step_matrix_mV_per_nA: np.ndarray,
        site_nodes: list[int] | tuple[int, ...] | np.ndarray,
        soma_node: int,
        *,
        rest_mV: float = -70.0,
        reversal_mV: float = 0.0,
        gamma_per_mV: float = 0.078,
        mg_divisor: float = 3.57,
        newton_max_iterations: int = 30,
        newton_tolerance_mV: float = 1e-10,
        max_backtrack: int = 8,
        min_alpha: float = 1.0 / 256.0,
    ) -> None:
        P = np.asarray(passive_step_matrix, dtype=float)
        X = np.asarray(input_step_matrix_mV_per_nA, dtype=float)
        if P.ndim != 2 or P.shape[0] != P.shape[1] or P.shape[0] == 0:
            raise ValueError("passive_step_matrix must be non-empty and square")
        n = P.shape[0]

        sites = np.asarray(site_nodes, dtype=int)
        if sites.ndim != 1 or len(sites) == 0:
            raise ValueError("site_nodes must be a non-empty 1-D sequence")
        if X.shape != (n, len(sites)):
            raise ValueError(
                f"input_step_matrix must have shape {(n, len(sites))}"
            )
        if np.any(sites < 0) or np.any(sites >= n):
            raise IndexError("site_nodes contain an invalid compartment index")
        if len(np.unique(sites)) != len(sites):
            raise ValueError("site_nodes must be unique")

        soma = int(soma_node)
        if soma < 0 or soma >= n:
            raise IndexError("soma_node is outside the compiled state")

        if not np.all(np.isfinite(P)) or not np.all(np.isfinite(X)):
            raise FloatingPointError("compiled operator contains non-finite values")

        self.P = P
        self.X_mV_per_nA = X
        self.site_nodes = sites
        self.soma_node = soma
        self.nstate = n
        self.nsite = len(sites)
        self.R_mV_per_nA = X[sites, :]

        self.rest_mV = float(rest_mV)
        self.reversal_mV = float(reversal_mV)
        self.gamma_per_mV = float(gamma_per_mV)
        self.mg_divisor = float(mg_divisor)

        self.newton_max_iterations = int(newton_max_iterations)
        self.newton_tolerance_mV = float(newton_tolerance_mV)
        self.max_backtrack = int(max_backtrack)
        self.min_alpha = float(min_alpha)
        if self.newton_max_iterations < 1:
            raise ValueError("newton_max_iterations must be positive")
        if self.max_backtrack < 0:
            raise ValueError("max_backtrack must be non-negative")
        if not (0.0 < self.min_alpha <= 1.0):
            raise ValueError("min_alpha must lie in (0, 1]")

    @classmethod
    def from_dense_passive_graph(
        cls,
        G_uS: np.ndarray,
        C_nF: np.ndarray,
        site_nodes: list[int] | tuple[int, ...] | np.ndarray,
        soma_node: int,
        *,
        dt_ms: float = 0.025,
        **kwargs,
    ) -> "CausalGraphCircuit":
        P, X = compile_dense_passive_graph(
            G_uS,
            C_nF,
            site_nodes,
            dt_ms=dt_ms,
        )
        return cls(P, X, site_nodes, soma_node, **kwargs)

    def current_law(
        self,
        voltage_mV: np.ndarray,
        g_ampa_uS: np.ndarray,
        g_nmda_raw_uS: np.ndarray,
    ) -> np.ndarray:
        return inward_synaptic_current_nA(
            voltage_mV,
            g_ampa_uS,
            g_nmda_raw_uS,
            reversal_mV=self.reversal_mV,
            gamma_per_mV=self.gamma_per_mV,
            mg_divisor=self.mg_divisor,
        )

    def current_derivative(
        self,
        voltage_mV: np.ndarray,
        g_ampa_uS: np.ndarray,
        g_nmda_raw_uS: np.ndarray,
    ) -> np.ndarray:
        return current_derivative_nA_per_mV(
            voltage_mV,
            g_ampa_uS,
            g_nmda_raw_uS,
            reversal_mV=self.reversal_mV,
            gamma_per_mV=self.gamma_per_mV,
            mg_divisor=self.mg_divisor,
        )

    def residual_and_jacobian(
        self,
        z_mV: np.ndarray,
        passive_sites_mV: np.ndarray,
        g_ampa_uS: np.ndarray,
        g_nmda_raw_uS: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        z = np.asarray(z_mV, dtype=float)
        passive = np.asarray(passive_sites_mV, dtype=float)
        absolute_v = self.rest_mV + z
        current = self.current_law(
            absolute_v,
            g_ampa_uS,
            g_nmda_raw_uS,
        )
        F = z - passive - self.R_mV_per_nA @ current
        dJ = self.current_derivative(
            absolute_v,
            g_ampa_uS,
            g_nmda_raw_uS,
        )
        jac = np.eye(self.nsite) - self.R_mV_per_nA @ np.diag(dJ)
        return F, jac, current

    def solve_site_step(
        self,
        passive_sites_mV: np.ndarray,
        g_ampa_uS: np.ndarray,
        g_nmda_raw_uS: np.ndarray,
        initial_site_depolarization_mV: np.ndarray,
    ) -> StepSolve:
        passive = np.asarray(passive_sites_mV, dtype=float)
        ga = np.asarray(g_ampa_uS, dtype=float)
        gn = np.asarray(g_nmda_raw_uS, dtype=float)
        z = np.asarray(initial_site_depolarization_mV, dtype=float).copy()
        expected = (self.nsite,)
        for name, value in (
            ("passive_sites_mV", passive),
            ("g_ampa_uS", ga),
            ("g_nmda_raw_uS", gn),
            ("initial_site_depolarization_mV", z),
        ):
            if value.shape != expected:
                raise ValueError(f"{name} must have shape {expected}")

        line_search_failures = 0
        for iteration in range(1, self.newton_max_iterations + 1):
            F, jac, current = self.residual_and_jacobian(
                z, passive, ga, gn
            )
            norm0 = float(np.max(np.abs(F)))
            if norm0 <= self.newton_tolerance_mV:
                return StepSolve(
                    site_depolarization_mV=z,
                    current_nA=current,
                    converged=True,
                    iterations=iteration - 1,
                    residual_inf_mV=norm0,
                    line_search_failures=line_search_failures,
                )

            try:
                delta = np.linalg.solve(jac, F)
            except np.linalg.LinAlgError:
                return StepSolve(
                    site_depolarization_mV=z,
                    current_nA=current,
                    converged=False,
                    iterations=iteration,
                    residual_inf_mV=norm0,
                    line_search_failures=line_search_failures + 1,
                    failure="singular_newton_jacobian",
                )

            accepted = False
            alpha = 1.0
            for _ in range(self.max_backtrack + 1):
                candidate = z - alpha * delta
                Fc, _, _ = self.residual_and_jacobian(
                    candidate, passive, ga, gn
                )
                if float(np.max(np.abs(Fc))) < norm0:
                    z = candidate
                    accepted = True
                    break
                alpha *= 0.5
                if alpha < self.min_alpha - 1e-15:
                    break

            if not accepted:
                line_search_failures += 1
                return StepSolve(
                    site_depolarization_mV=z,
                    current_nA=current,
                    converged=False,
                    iterations=iteration,
                    residual_inf_mV=norm0,
                    line_search_failures=line_search_failures,
                    failure="backtracking_no_decrease",
                )

        F, _, current = self.residual_and_jacobian(
            z, passive, ga, gn
        )
        return StepSolve(
            site_depolarization_mV=z,
            current_nA=current,
            converged=False,
            iterations=self.newton_max_iterations,
            residual_inf_mV=float(np.max(np.abs(F))),
            line_search_failures=line_search_failures,
            failure="newton_iteration_limit",
        )

    def run(
        self,
        g_ampa_uS: np.ndarray,
        g_nmda_raw_uS: np.ndarray,
        *,
        initial_state_depolarization_mV: np.ndarray | None = None,
    ) -> CausalSolveResult:
        ga = np.asarray(g_ampa_uS, dtype=float)
        gn = np.asarray(g_nmda_raw_uS, dtype=float)
        if ga.ndim != 2 or ga.shape[0] != self.nsite:
            raise ValueError(
                f"g_ampa_uS must have shape ({self.nsite}, time)"
            )
        if gn.shape != ga.shape:
            raise ValueError("g_nmda_raw_uS must match g_ampa_uS")
        if not np.all(np.isfinite(ga)) or not np.all(np.isfinite(gn)):
            raise FloatingPointError("conductance arrays contain non-finite values")
        if np.any(ga < 0) or np.any(gn < 0):
            raise ValueError("conductances must be non-negative")

        ntime = ga.shape[1]
        if initial_state_depolarization_mV is None:
            state = np.zeros(self.nstate, dtype=float)
        else:
            state = np.asarray(
                initial_state_depolarization_mV, dtype=float
            ).copy()
            if state.shape != (self.nstate,):
                raise ValueError(
                    f"initial state must have shape {(self.nstate,)}"
                )

        previous_z = state[self.site_nodes].copy()
        pred_current = np.zeros((self.nsite, ntime), dtype=float)
        pred_local_abs = np.zeros((self.nsite, ntime), dtype=float)
        pred_soma = np.zeros(ntime, dtype=float)

        all_converged = True
        max_iterations = 0
        max_residual = 0.0
        line_search_failures = 0
        max_consistency = 0.0
        first_failure = None

        for ti in range(ntime):
            passive = self.P @ state
            passive_sites = passive[self.site_nodes]
            solved = self.solve_site_step(
                passive_sites,
                ga[:, ti],
                gn[:, ti],
                previous_z,
            )
            max_iterations = max(max_iterations, solved.iterations)
            max_residual = max(max_residual, solved.residual_inf_mV)
            line_search_failures += solved.line_search_failures
            if not solved.converged:
                all_converged = False
                if first_failure is None:
                    first_failure = ti

            z = solved.site_depolarization_mV
            current = self.current_law(
                self.rest_mV + z,
                ga[:, ti],
                gn[:, ti],
            )
            state = passive + self.X_mV_per_nA @ current
            consistency = float(
                np.max(np.abs(state[self.site_nodes] - z))
            )
            max_consistency = max(max_consistency, consistency)

            pred_current[:, ti] = current
            pred_local_abs[:, ti] = self.rest_mV + state[self.site_nodes]
            pred_soma[ti] = state[self.soma_node]
            previous_z = state[self.site_nodes].copy()

        return CausalSolveResult(
            current_nA=pred_current,
            local_voltage_mV=pred_local_abs,
            soma_depolarization_mV=pred_soma,
            final_state_depolarization_mV=state.copy(),
            all_steps_converged=all_converged,
            max_newton_iterations=max_iterations,
            max_residual_inf_mV=max_residual,
            total_line_search_failures=line_search_failures,
            max_site_consistency_mV=max_consistency,
            first_failure_time_index=first_failure,
        )
