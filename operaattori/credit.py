"""Gate-4 delayed structural credit.

This module is intentionally separate from the Gate-1 order-memory substrate.
The question is narrower: can a delayed scalar soma consequence act through a
lingering *local* eligibility field and make the final material operator more
selective?

The target is deliberately easy. An explicit two-weight digital model solves
it immediately. The purpose is to test credit/addressability in grown
material, not benchmark difficulty.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CreditConfig:
    n: int = 32
    dt: float = 0.02

    base_conductance: float = 0.20
    morph_conductance: float = 8.0
    diffusion: float = 0.80
    fast_decay: float = 0.12
    input_gain: float = 1.5

    eligibility_decay: float = 0.08
    eligibility_gain: float = 1.0

    drive_steps: int = 50
    observe_steps: int = 100
    credit_delay_steps: int = 30
    recovery_steps: int = 80
    trials: int = 120

    learning_rate: float = 5.0
    trial_morph_decay: float = 0.001

    target_a: float = 0.008
    target_b: float = 0.0002

    eval_steps: int = 800
    eval_pulse_steps: int = 60


def _field_step(
    fast: np.ndarray,
    eligibility: np.ndarray,
    morphology: np.ndarray,
    stimulus: tuple[float, float],
    config: CreditConfig,
) -> tuple[np.ndarray, np.ndarray]:
    node_g = (
        config.base_conductance
        + config.morph_conductance * morphology
    )
    edge_g = 0.5 * (node_g[:-1] + node_g[1:])
    flux = edge_g * (fast[:-1] - fast[1:])

    divergence = np.zeros_like(fast)
    divergence[:-1] -= flux
    divergence[1:] += flux

    du = config.diffusion * divergence - config.fast_decay * fast
    du[2] += config.input_gain * stimulus[0]
    du[-3] += config.input_gain * stimulus[1]
    fast = fast + config.dt * du

    local_transport = np.zeros(config.n, dtype=float)
    local_transport[:-1] += np.abs(flux)
    local_transport[1:] += np.abs(flux)
    local_transport *= 0.5
    eligibility = eligibility + config.dt * (
        config.eligibility_gain * local_transport
        - config.eligibility_decay * eligibility
    )
    return fast, eligibility


def evaluate_ports(
    morphology: np.ndarray,
    config: CreditConfig | None = None,
) -> tuple[float, float]:
    c = config or CreditConfig()
    out: list[float] = []
    soma = c.n // 2
    for port in (0, 1):
        fast = np.zeros(c.n, dtype=float)
        eligibility = np.zeros(c.n, dtype=float)
        peak = 0.0
        for t in range(c.eval_steps):
            if t < c.eval_pulse_steps:
                stimulus = (1.0, 0.0) if port == 0 else (0.0, 1.0)
            else:
                stimulus = (0.0, 0.0)
            fast, eligibility = _field_step(
                fast, eligibility, morphology, stimulus, c
            )
            peak = max(peak, float(fast[soma]))
        out.append(peak)
    return out[0], out[1]


def selectivity(a_response: float, b_response: float) -> float:
    return (
        (a_response - b_response)
        / (a_response + b_response + 1e-12)
    )


def train(
    arm: str,
    seed: int,
    config: CreditConfig | None = None,
) -> dict[str, object]:
    """Train one continuous substrate.

    Arms:
      causal               correct delayed scalar × local eligibility
      credit_shuffle       same credit magnitude, randomized sign
      eligibility_shuffle  correct credit × spatially permuted eligibility
      no_credit            no structural consequence
    """
    c = config or CreditConfig()
    if arm not in {
        "causal",
        "credit_shuffle",
        "eligibility_shuffle",
        "no_credit",
    }:
        raise ValueError(arm)

    init_rng = np.random.default_rng(seed)
    credit_rng = np.random.default_rng(50_000 + seed)
    spatial_rng = np.random.default_rng(90_000 + seed)

    morphology = np.clip(
        0.05 + 0.005 * init_rng.normal(size=c.n),
        0.0,
        1.0,
    )
    fast = np.zeros(c.n, dtype=float)
    eligibility = np.zeros(c.n, dtype=float)

    trial_ports = np.asarray([0, 1] * (c.trials // 2), dtype=int)
    np.random.default_rng(999 + seed).shuffle(trial_ports)

    consequences: list[float] = []

    for port in trial_ports:
        peak = 0.0
        soma = c.n // 2

        # Drive.
        for _ in range(c.drive_steps):
            stimulus = (1.0, 0.0) if port == 0 else (0.0, 1.0)
            fast, eligibility = _field_step(
                fast, eligibility, morphology, stimulus, c
            )
            peak = max(peak, float(fast[soma]))

        # Observe the consequence during silence.
        for _ in range(c.observe_steps):
            fast, eligibility = _field_step(
                fast,
                eligibility,
                morphology,
                (0.0, 0.0),
                c,
            )
            peak = max(peak, float(fast[soma]))

        target = c.target_a if port == 0 else c.target_b
        consequence = target - peak
        consequences.append(consequence)

        # Consequence is withheld while local eligibility continues to decay.
        for _ in range(c.credit_delay_steps):
            fast, eligibility = _field_step(
                fast,
                eligibility,
                morphology,
                (0.0, 0.0),
                c,
            )

        if arm == "causal":
            scalar = consequence
            addressed = eligibility
        elif arm == "credit_shuffle":
            scalar = abs(consequence) * (
                1.0 if credit_rng.random() < 0.5 else -1.0
            )
            addressed = eligibility
        elif arm == "eligibility_shuffle":
            scalar = consequence
            addressed = eligibility[
                spatial_rng.permutation(c.n)
            ]
        else:
            scalar = 0.0
            addressed = eligibility

        morphology = np.clip(
            morphology
            + c.learning_rate * scalar * addressed
            - c.trial_morph_decay * morphology,
            0.0,
            1.0,
        )

        # The machine is not reset. It simply receives more silence.
        for _ in range(c.recovery_steps):
            fast, eligibility = _field_step(
                fast,
                eligibility,
                morphology,
                (0.0, 0.0),
                c,
            )

    a_response, b_response = evaluate_ports(morphology, c)
    return {
        "morphology": morphology,
        "a_response": a_response,
        "b_response": b_response,
        "selectivity": selectivity(a_response, b_response),
        "mean_abs_consequence": float(
            np.mean(np.abs(consequences))
        ),
    }
