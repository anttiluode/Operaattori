"""Portable Green-matrix nonlinear circuit runtime.

This module contains no NEURON or FCI dependency.  It implements the reduced
operator earned by Operaattori's held-out audits:

    local voltage = baseline + G * current
    current       = conductance_law(local voltage)
    soma voltage  = T * current

The caller supplies measured linear kernels and conductance templates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _next_pow_two(n: int) -> int:
    return 1 << int(np.ceil(np.log2(max(1, int(n)))))


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
    """Return positive-inward AMPA+NMDA current in nA.

    uS * mV = nA.  The raw NMDA conductance excludes magnesium block.
    """
    v = np.asarray(voltage_mV, dtype=float)
    ga = np.asarray(g_ampa_uS, dtype=float)
    gn = np.asarray(g_nmda_raw_uS, dtype=float)
    block = magnesium_block(
        v,
        gamma_per_mV=gamma_per_mV,
        mg_divisor=mg_divisor,
    )
    return (ga + gn * block) * (float(reversal_mV) - v)


def shift_template(
    template: np.ndarray,
    delay_samples: int,
) -> np.ndarray:
    """Causally shift a finite conductance template without wraparound."""
    x = np.asarray(template, dtype=float)
    k = int(delay_samples)
    if k < 0:
        raise ValueError("delay_samples must be non-negative")
    out = np.zeros_like(x)
    if k < len(x):
        out[k:] = x[: len(x) - k]
    return out


def timed_conductance_matrix(
    template: np.ndarray,
    delay_samples: list[int] | tuple[int, ...] | np.ndarray,
) -> np.ndarray:
    """Build one shifted copy of a conductance template per site."""
    return np.stack(
        [shift_template(template, int(k)) for k in delay_samples],
        axis=0,
    )


@dataclass(frozen=True)
class SolveResult:
    current_nA: np.ndarray
    local_voltage_mV: np.ndarray
    soma_depolarization_mV: np.ndarray
    converged: bool
    iterations: int
    final_relative_current_update: float


class GreenCircuit:
    """Nonlinear point-conductance circuit embedded in linear Green kernels.

    Parameters
    ----------
    local_kernel_mV_per_nA_sample:
        Array [target_site, source_site, time].  Each kernel is the local
        voltage response to a +1 nA one-sample current at the source site.
    soma_kernel_mV_per_nA_sample:
        Array [source_site, time] with the corresponding soma response.
    baseline_local_mV:
        Array [site, time] giving the no-input local-voltage trajectory.
    """

    def __init__(
        self,
        local_kernel_mV_per_nA_sample: np.ndarray,
        soma_kernel_mV_per_nA_sample: np.ndarray,
        baseline_local_mV: np.ndarray,
        *,
        gamma_per_mV: float = 0.078,
        mg_divisor: float = 3.57,
        reversal_mV: float = 0.0,
        damping: float = 0.5,
        relative_tolerance: float = 1e-8,
        max_iterations: int = 200,
    ) -> None:
        local_h = np.asarray(
            local_kernel_mV_per_nA_sample,
            dtype=float,
        )
        soma_h = np.asarray(
            soma_kernel_mV_per_nA_sample,
            dtype=float,
        )
        baseline = np.asarray(
            baseline_local_mV,
            dtype=float,
        )

        if local_h.ndim != 3:
            raise ValueError("local kernel must have shape [site, site, time]")
        if soma_h.ndim != 2:
            raise ValueError("soma kernel must have shape [site, time]")
        if baseline.ndim != 2:
            raise ValueError("baseline must have shape [site, time]")
        nsite, nsource, ntime = local_h.shape
        if nsite != nsource:
            raise ValueError("local kernel site axes must be square")
        if soma_h.shape != (nsite, ntime):
            raise ValueError("soma kernel shape does not match local kernel")
        if baseline.shape != (nsite, ntime):
            raise ValueError("baseline shape does not match kernel")
        if not (0.0 < float(damping) <= 1.0):
            raise ValueError("damping must lie in (0, 1]")
        if int(max_iterations) < 1:
            raise ValueError("max_iterations must be positive")

        self.local_h = local_h
        self.soma_h = soma_h
        self.baseline = baseline
        self.nsite = nsite
        self.ntime = ntime
        self.gamma = float(gamma_per_mV)
        self.mg_divisor = float(mg_divisor)
        self.reversal = float(reversal_mV)
        self.damping = float(damping)
        self.tol = float(relative_tolerance)
        self.max_iterations = int(max_iterations)

        self.nfft = _next_pow_two(2 * ntime - 1)
        self._local_f = np.fft.rfft(
            local_h,
            n=self.nfft,
            axis=-1,
        )
        self._soma_f = np.fft.rfft(
            soma_h,
            n=self.nfft,
            axis=-1,
        )

    @staticmethod
    def _rms(x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        return float(np.sqrt(np.mean(x * x)))

    def local_transport(self, current_nA: np.ndarray) -> np.ndarray:
        current = np.asarray(current_nA, dtype=float)
        if current.shape != (self.nsite, self.ntime):
            raise ValueError("current shape does not match circuit")
        jf = np.fft.rfft(
            current,
            n=self.nfft,
            axis=-1,
        )
        vf = np.einsum(
            "ijw,jw->iw",
            self._local_f,
            jf,
            optimize=True,
        )
        return np.fft.irfft(
            vf,
            n=self.nfft,
            axis=-1,
        )[:, : self.ntime]

    def soma_transport(self, current_nA: np.ndarray) -> np.ndarray:
        current = np.asarray(current_nA, dtype=float)
        if current.shape != (self.nsite, self.ntime):
            raise ValueError("current shape does not match circuit")
        jf = np.fft.rfft(
            current,
            n=self.nfft,
            axis=-1,
        )
        vf = np.sum(self._soma_f * jf, axis=0)
        return np.fft.irfft(
            vf,
            n=self.nfft,
        )[: self.ntime]

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
            reversal_mV=self.reversal,
            gamma_per_mV=self.gamma,
            mg_divisor=self.mg_divisor,
        )

    def solve(
        self,
        g_ampa_uS: np.ndarray,
        g_nmda_raw_uS: np.ndarray,
    ) -> SolveResult:
        ga = np.asarray(g_ampa_uS, dtype=float)
        gn = np.asarray(g_nmda_raw_uS, dtype=float)
        expected = (self.nsite, self.ntime)
        if ga.shape != expected or gn.shape != expected:
            raise ValueError(
                f"conductance arrays must have shape {expected}"
            )
        if not np.all(np.isfinite(ga)) or not np.all(np.isfinite(gn)):
            raise FloatingPointError("non-finite conductance")

        current = self.current_law(
            self.baseline,
            ga,
            gn,
        )
        converged = False
        final_error = float("inf")
        iterations = 0

        for iteration in range(1, self.max_iterations + 1):
            voltage = self.baseline + self.local_transport(current)
            target = self.current_law(voltage, ga, gn)
            updated = (
                (1.0 - self.damping) * current
                + self.damping * target
            )
            final_error = self._rms(updated - current) / (
                self._rms(updated) + 1e-30
            )
            current = updated
            iterations = iteration
            if final_error <= self.tol:
                converged = True
                break

        voltage = self.baseline + self.local_transport(current)
        soma = self.soma_transport(current)
        return SolveResult(
            current_nA=current,
            local_voltage_mV=voltage,
            soma_depolarization_mV=soma,
            converged=converged,
            iterations=iterations,
            final_relative_current_update=float(final_error),
        )
