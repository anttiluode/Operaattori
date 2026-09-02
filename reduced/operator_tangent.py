"""Exact first-order tangents through the compiled causal graph circuit.

This module differentiates the object Operaattori actually earned:

    (G, C) -> (P, X) -> local implicit NMDA -> state -> soma

No autodiff package is required.  The passive compiler derivative follows from
A^-1 differentiation and the nonlinear-site derivative follows from the same
implicit Jacobian already used by the causal Newton solve.

The API supports several scalar geometry/control parameters at once.  A caller
supplies dG/dtheta and dC/dtheta for each parameter; conductance programs are
held fixed.  This is deliberately a tangent of the compiled circuit, not a
claim that any particular dendritic shape objective is biologically preferred.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from reduced.causal_graph_circuit import CausalGraphCircuit


@dataclass(frozen=True)
class CompiledOperatorTangent:
    passive_step_matrix: np.ndarray
    input_step_matrix_mV_per_nA: np.ndarray
    dpassive_step_dtheta: np.ndarray
    dinput_step_dtheta_mV_per_nA: np.ndarray


@dataclass(frozen=True)
class CausalTangentResult:
    soma_depolarization_mV: np.ndarray
    soma_tangent_mV_per_unit: np.ndarray
    final_state_depolarization_mV: np.ndarray
    final_state_tangent_mV_per_unit: np.ndarray
    all_steps_converged: bool
    max_newton_iterations: int
    max_site_consistency_mV: float
    max_site_tangent_consistency_mV_per_unit: float


def _as_parameter_stack(value: np.ndarray, trailing_shape: tuple[int, ...], name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape == trailing_shape:
        arr = arr[None, ...]
    if arr.ndim != len(trailing_shape) + 1 or arr.shape[1:] != trailing_shape:
        raise ValueError(f"{name} must have shape {trailing_shape} or (parameter, {', '.join(map(str, trailing_shape))})")
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return arr


def compile_dense_passive_graph_tangent(
    G_uS: np.ndarray,
    C_nF: np.ndarray,
    dG_dtheta_uS: np.ndarray,
    dC_dtheta_nF: np.ndarray,
    site_nodes: list[int] | tuple[int, ...] | np.ndarray,
    *,
    dt_ms: float = 0.025,
) -> CompiledOperatorTangent:
    """Compile P, X and exact first derivatives for one or more parameters.

    For A = G + D, D = diag(C/dt), P = A^-1 D, X = A^-1 B:

        dP = A^-1 (dD - dA P)
        dX = -A^-1 dA X

    The helper is dense and intended for small/reference graphs.  Large
    morphologies can use the same equations with sparse linear solves.
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
    if np.any(C <= 0):
        raise ValueError("C_nF must be strictly positive")

    dG = _as_parameter_stack(dG_dtheta_uS, (n, n), "dG_dtheta_uS")
    dC = _as_parameter_stack(dC_dtheta_nF, (n,), "dC_dtheta_nF")
    if dG.shape[0] != dC.shape[0]:
        raise ValueError("dG_dtheta_uS and dC_dtheta_nF must have the same parameter count")

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
    P = solved[:, :n]
    X = solved[:, n:]

    nparam = dG.shape[0]
    dP = np.empty((nparam, n, n), dtype=float)
    dX = np.empty((nparam, n, len(sites)), dtype=float)
    for k in range(nparam):
        dd = dC[k] / float(dt_ms)
        dD = np.diag(dd)
        dA = dG[k] + dD
        tangent_rhs = np.concatenate(
            (dD - dA @ P, -(dA @ X)),
            axis=1,
        )
        tangent = np.linalg.solve(A, tangent_rhs)
        dP[k] = tangent[:, :n]
        dX[k] = tangent[:, n:]

    return CompiledOperatorTangent(P, X, dP, dX)


def run_operator_tangent(
    circuit: CausalGraphCircuit,
    g_ampa_uS: np.ndarray,
    g_nmda_raw_uS: np.ndarray,
    dpassive_step_dtheta: np.ndarray,
    dinput_step_dtheta_mV_per_nA: np.ndarray,
) -> CausalTangentResult:
    """Run the causal circuit and propagate exact operator tangents."""
    ga = np.asarray(g_ampa_uS, dtype=float)
    gn = np.asarray(g_nmda_raw_uS, dtype=float)
    if ga.ndim != 2 or ga.shape[0] != circuit.nsite:
        raise ValueError(f"g_ampa_uS must have shape ({circuit.nsite}, time)")
    if gn.shape != ga.shape:
        raise ValueError("g_nmda_raw_uS must match g_ampa_uS")
    if np.any(ga < 0) or np.any(gn < 0):
        raise ValueError("conductances must be non-negative")

    dP = _as_parameter_stack(
        dpassive_step_dtheta,
        (circuit.nstate, circuit.nstate),
        "dpassive_step_dtheta",
    )
    dX = _as_parameter_stack(
        dinput_step_dtheta_mV_per_nA,
        (circuit.nstate, circuit.nsite),
        "dinput_step_dtheta_mV_per_nA",
    )
    if dP.shape[0] != dX.shape[0]:
        raise ValueError("dP and dX must have the same parameter count")

    nparam = dP.shape[0]
    ntime = ga.shape[1]
    state = np.zeros(circuit.nstate, dtype=float)
    dstate = np.zeros((nparam, circuit.nstate), dtype=float)
    previous_z = np.zeros(circuit.nsite, dtype=float)

    soma = np.zeros(ntime, dtype=float)
    soma_tangent = np.zeros((nparam, ntime), dtype=float)
    all_converged = True
    max_iterations = 0
    max_consistency = 0.0
    max_tangent_consistency = 0.0

    R = circuit.R_mV_per_nA
    dR = dX[:, circuit.site_nodes, :]

    for ti in range(ntime):
        old_state = state
        old_dstate = dstate

        passive = circuit.P @ old_state
        dpassive = np.empty_like(old_dstate)
        for k in range(nparam):
            dpassive[k] = dP[k] @ old_state + circuit.P @ old_dstate[k]

        passive_sites = passive[circuit.site_nodes]
        solved = circuit.solve_site_step(
            passive_sites,
            ga[:, ti],
            gn[:, ti],
            previous_z,
        )
        all_converged = all_converged and solved.converged
        max_iterations = max(max_iterations, solved.iterations)

        z = solved.site_depolarization_mV
        absolute_v = circuit.rest_mV + z
        current = circuit.current_law(absolute_v, ga[:, ti], gn[:, ti])
        dJdV = circuit.current_derivative(absolute_v, ga[:, ti], gn[:, ti])
        K = np.eye(circuit.nsite) - R @ np.diag(dJdV)

        dz = np.empty((nparam, circuit.nsite), dtype=float)
        dcurrent = np.empty((nparam, circuit.nsite), dtype=float)
        for k in range(nparam):
            rhs = dpassive[k, circuit.site_nodes] + dR[k] @ current
            dz[k] = np.linalg.solve(K, rhs)
            dcurrent[k] = dJdV * dz[k]

        state = passive + circuit.X_mV_per_nA @ current
        dstate = np.empty_like(old_dstate)
        for k in range(nparam):
            dstate[k] = (
                dpassive[k]
                + dX[k] @ current
                + circuit.X_mV_per_nA @ dcurrent[k]
            )

        max_consistency = max(
            max_consistency,
            float(np.max(np.abs(state[circuit.site_nodes] - z))),
        )
        max_tangent_consistency = max(
            max_tangent_consistency,
            float(np.max(np.abs(dstate[:, circuit.site_nodes] - dz))),
        )
        soma[ti] = state[circuit.soma_node]
        soma_tangent[:, ti] = dstate[:, circuit.soma_node]
        previous_z = state[circuit.site_nodes].copy()

    return CausalTangentResult(
        soma_depolarization_mV=soma,
        soma_tangent_mV_per_unit=soma_tangent,
        final_state_depolarization_mV=state.copy(),
        final_state_tangent_mV_per_unit=dstate.copy(),
        all_steps_converged=all_converged,
        max_newton_iterations=max_iterations,
        max_site_consistency_mV=max_consistency,
        max_site_tangent_consistency_mV_per_unit=max_tangent_consistency,
    )
