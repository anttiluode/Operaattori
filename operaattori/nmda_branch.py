"""Reduced, source-grounded NMDA nonlinearity on the real passive tree.

This is a quasi-static peak-conductance mechanism audit, not a reproduction of
Aizenbud et al. Figure 4 or their full NEURON model. The passive surrounding
tree is exact at DC; selected synaptic sites close a self-consistent local
voltage feedback loop through the Jahr-Stevens magnesium gate used by the
paper and released FCI code.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SynapseCondition:
    name: str
    g_ampa_uS: float
    g_nmda_uS: float
    gamma_per_mV: float
    frozen_block: bool = False
    fixed_current: bool = False


HUMAN = SynapseCondition("human", 0.00088, 0.00131, 0.078)
HYBRID_B = SynapseCondition("hybrid_b", 0.00088, 0.00131, 0.062)
HUMAN_FROZEN_BLOCK = SynapseCondition(
    "human_frozen_block", 0.00088, 0.00131, 0.078, frozen_block=True
)
HUMAN_AMPA_ONLY = SynapseCondition("human_ampa_only", 0.00088, 0.0, 0.078)
HUMAN_LINEAR_CURRENT = SynapseCondition(
    "human_linear_current", 0.00088, 0.00131, 0.078, fixed_current=True
)


def nmda_block(v_mV: np.ndarray | float, gamma_per_mV: float) -> np.ndarray:
    """Jahr-Stevens block used by Aizenbud et al.; Mg=1 mM, n=1/3.57 mM^-1."""
    v = np.asarray(v_mV, dtype=float)
    return 1.0 / (1.0 + np.exp(-float(gamma_per_mV) * v) * (1.0 / 3.57))


def synaptic_current_and_derivative_nA(
    v_mV: np.ndarray,
    multiplicity: np.ndarray,
    condition: SynapseCondition,
    *,
    v_rest_mV: float = -70.0,
    e_syn_mV: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Inward current and dI/dV at selected sites, in nA and nA/mV."""
    v = np.asarray(v_mV, dtype=float)
    mult = np.asarray(multiplicity, dtype=float)
    if condition.fixed_current:
        b0 = float(nmda_block(v_rest_mV, condition.gamma_per_mV))
        geff = condition.g_ampa_uS + condition.g_nmda_uS * b0
        current = mult * geff * (e_syn_mV - v_rest_mV)
        return current, np.zeros_like(current)

    if condition.frozen_block:
        b = np.full_like(v, float(nmda_block(v_rest_mV, condition.gamma_per_mV)))
        db = np.zeros_like(v)
    else:
        b = nmda_block(v, condition.gamma_per_mV)
        db = condition.gamma_per_mV * b * (1.0 - b)

    geff = condition.g_ampa_uS + condition.g_nmda_uS * b
    current = mult * geff * (e_syn_mV - v)
    deriv = mult * (
        condition.g_nmda_uS * db * (e_syn_mV - v) - geff
    )
    return current, deriv


def solve_equilibrium(
    z_mohm: np.ndarray,
    multiplicity: np.ndarray,
    condition: SynapseCondition,
    *,
    v_rest_mV: float = -70.0,
    max_iter: int = 80,
    tol_mV: float = 1e-9,
) -> dict:
    """Solve v = v_rest + Z I_syn(v) by damped Newton iteration."""
    Z = np.asarray(z_mohm, dtype=float)
    mult = np.asarray(multiplicity, dtype=float)
    m = len(mult)
    if Z.shape != (m, m):
        raise ValueError("Z and multiplicity size mismatch")

    dep = np.zeros(m, dtype=float)
    eye = np.eye(m)
    converged = False
    residual = np.inf

    for it in range(max_iter):
        absolute_v = v_rest_mV + dep
        current, deriv = synaptic_current_and_derivative_nA(
            absolute_v, mult, condition, v_rest_mV=v_rest_mV
        )
        F = dep - Z @ current
        residual = float(np.max(np.abs(F))) if m else 0.0
        if residual < tol_mV:
            converged = True
            break
        J = eye - Z * deriv[np.newaxis, :]
        try:
            step = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            step = -0.1 * F

        # Backtracking keeps Newton on the physically continuous low-voltage
        # branch until that branch itself becomes unstable.
        accepted = False
        scale = 1.0
        for _ in range(12):
            trial = np.clip(dep + scale * step, -20.0, 69.9)
            tv = v_rest_mV + trial
            ti, _ = synaptic_current_and_derivative_nA(
                tv, mult, condition, v_rest_mV=v_rest_mV
            )
            tres = float(np.max(np.abs(trial - Z @ ti)))
            if tres < residual:
                dep = trial
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            dep = np.clip(dep - 0.05 * F, -20.0, 69.9)

    absolute_v = v_rest_mV + dep
    current, deriv = synaptic_current_and_derivative_nA(
        absolute_v, mult, condition, v_rest_mV=v_rest_mV
    )
    residual = float(np.max(np.abs(dep - Z @ current))) if m else 0.0
    converged = converged or residual < 1e-6
    return {
        "voltage_mV": absolute_v,
        "depolarization_mV": dep,
        "current_nA": current,
        "dI_dV_nA_per_mV": deriv,
        "converged": bool(converged),
        "residual_mV": residual,
    }


def effective_rank_rows(curves: np.ndarray) -> float:
    """Participation rank after normalizing each dose-response curve at dose 1."""
    X = np.asarray(curves, dtype=float)
    if X.ndim != 2 or X.shape[1] < 2:
        return 0.0
    denom = np.maximum(np.abs(X[:, [0]]), 1e-15)
    X = X / denom
    X = X - np.mean(X, axis=0, keepdims=True)
    if np.allclose(X, 0):
        return 0.0
    s = np.linalg.svd(X, compute_uv=False)
    p = s * s
    d = float(np.sum(p * p))
    return float((np.sum(p) ** 2) / d) if d > 0 else 0.0


def linear_regression_r2(X: np.ndarray, y: np.ndarray) -> float:
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    A = np.column_stack([np.ones(len(X)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
