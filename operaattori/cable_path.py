"""Passive cable transfer matrices on a branching matrix scaffold.

This module is deliberately small and classical.  Each cylindrical dendritic
segment is represented as a frequency-domain two-port (ABCD) matrix for the
passive cable equation:

    [V_prox]   [A B] [V_dist]
    [I_prox] = [C D] [I_dist]

For a uniform segment

    gamma = sqrt(r_a * y_m)
    Z0    = sqrt(r_a / y_m)

    M(l) = [[cosh(gamma l), Z0 sinh(gamma l)],
            [sinh(gamma l)/Z0, cosh(gamma l)]]

where r_a is axial resistance per unit length and y_m is membrane admittance
per unit length.

At fixed frequency, segments with the same radius share the same generator and
M(l1) M(l2) = M(l1+l2): order is irrelevant.  Different radii change both
gamma and Z0, and the corresponding local transport matrices generally do not
commute.  This makes a clean Gate-11 assay for whether *spatial attachment and
order* of local operators can matter on the real neuron.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class PassiveCableParams:
    # Aizenbud et al. 2026 detailed models use these passive values.
    rm_ohm_cm2: float = 20_000.0
    ra_ohm_cm: float = 150.0
    cm_uF_cm2: float = 1.0

    @property
    def cm_F_cm2(self) -> float:
        return self.cm_uF_cm2 * 1e-6


def cable_constants(
    radius_um: float,
    frequency_hz: float,
    params: PassiveCableParams | None = None,
) -> tuple[complex, complex]:
    """Return propagation constant gamma [1/cm] and characteristic Z0 [ohm]."""
    p = params or PassiveCableParams()
    a_cm = max(float(radius_um), 1e-6) * 1e-4
    omega = 2.0 * np.pi * float(frequency_hz)

    r_a = p.ra_ohm_cm / (np.pi * a_cm * a_cm)  # ohm / cm
    y_m = 2.0 * np.pi * a_cm * (
        (1.0 / p.rm_ohm_cm2) + 1j * omega * p.cm_F_cm2
    )  # siemens / cm

    gamma = complex(np.sqrt(r_a * y_m))
    z0 = complex(np.sqrt(r_a / y_m))
    return gamma, z0


def cable_abcd(
    length_um: float,
    radius_um: float,
    frequency_hz: float,
    params: PassiveCableParams | None = None,
) -> np.ndarray:
    """ABCD matrix mapping distal voltage/current to proximal voltage/current."""
    l_cm = max(float(length_um), 0.0) * 1e-4
    gamma, z0 = cable_constants(radius_um, frequency_hz, params)
    gl = gamma * l_cm
    ch = np.cosh(gl)
    sh = np.sinh(gl)
    return np.asarray(
        [[ch, z0 * sh], [sh / z0, ch]],
        dtype=np.complex128,
    )


def compose_abcd(matrices: Sequence[np.ndarray]) -> np.ndarray:
    """Compose matrices in proximal -> distal segment order."""
    out = np.eye(2, dtype=np.complex128)
    for M in matrices:
        out = out @ np.asarray(M, dtype=np.complex128)
    return out


def path_abcd(
    lengths_um: Sequence[float],
    radii_um: Sequence[float],
    frequency_hz: float,
    params: PassiveCableParams | None = None,
    order: Sequence[int] | None = None,
) -> np.ndarray:
    lengths = np.asarray(lengths_um, dtype=float)
    radii = np.asarray(radii_um, dtype=float)
    if lengths.shape != radii.shape:
        raise ValueError("lengths and radii must have the same shape")
    idx = np.arange(len(lengths)) if order is None else np.asarray(order, dtype=int)
    mats = [
        cable_abcd(lengths[i], radii[i], frequency_hz, params)
        for i in idx
    ]
    return compose_abcd(mats)


def sealed_input_impedance(M: np.ndarray) -> complex:
    """Input impedance when the distal end is sealed (I_dist = 0)."""
    M = np.asarray(M, dtype=np.complex128)
    C = M[1, 0]
    if abs(C) < 1e-300:
        return complex(np.inf)
    return complex(M[0, 0] / C)


def sealed_distal_voltage_gain(M: np.ndarray) -> complex:
    """V_dist / V_prox under sealed distal boundary."""
    A = complex(np.asarray(M, dtype=np.complex128)[0, 0])
    if abs(A) < 1e-300:
        return complex(np.inf)
    return 1.0 / A


def relative_complex_difference(a: complex, b: complex, eps: float = 1e-30) -> float:
    return float(abs(a - b) / (0.5 * (abs(a) + abs(b)) + eps))


def phase_difference(a: complex, b: complex) -> float:
    """Absolute wrapped phase difference in radians."""
    return float(abs(np.angle(a * np.conj(b))))


def area_preserving_radius(lengths_um: Sequence[float], radii_um: Sequence[float]) -> float:
    """Uniform radius preserving sum(2*pi*r*l) for the path."""
    l = np.asarray(lengths_um, dtype=float)
    r = np.asarray(radii_um, dtype=float)
    total = float(np.sum(l))
    if total <= 0:
        return float(np.mean(r))
    return float(np.sum(l * r) / total)


def commutator_action_score(
    A: np.ndarray,
    B: np.ndarray,
    *,
    test_states: Iterable[np.ndarray] | None = None,
) -> float:
    """Dimensionless, state-action measure of AB != BA.

    Raw ABCD matrix norms mix voltage/current units.  Instead, act on a few
    normalized [V, I] states after scaling current by a characteristic
    impedance inferred from A/B.  The score is zero iff the tested actions
    commute on those states.
    """
    A = np.asarray(A, dtype=np.complex128)
    B = np.asarray(B, dtype=np.complex128)

    # Use a robust scale from the B/C entries.  This only nondimensionalizes
    # the state for the diagnostic; functional assays below use physical Z/gain.
    z_candidates = []
    for M in (A, B):
        if abs(M[1, 0]) > 1e-30 and abs(M[0, 1]) > 1e-30:
            z_candidates.append(float(np.sqrt(abs(M[0, 1] / M[1, 0]))))
    z = float(np.median(z_candidates)) if z_candidates else 1.0
    z = max(z, 1e-12)

    P = np.asarray([[1.0, 0.0], [0.0, z]], dtype=np.complex128)
    Pinv = np.asarray([[1.0, 0.0], [0.0, 1.0 / z]], dtype=np.complex128)
    Ad = P @ A @ Pinv
    Bd = P @ B @ Pinv

    states = list(test_states) if test_states is not None else [
        np.asarray([1.0, 0.0], dtype=np.complex128),
        np.asarray([0.0, 1.0], dtype=np.complex128),
        np.asarray([1.0, 1.0j], dtype=np.complex128) / np.sqrt(2.0),
    ]
    vals = []
    for x in states:
        ab = Ad @ (Bd @ x)
        ba = Bd @ (Ad @ x)
        vals.append(float(np.linalg.norm(ab - ba) / (0.5 * (np.linalg.norm(ab) + np.linalg.norm(ba)) + 1e-30)))
    return float(np.mean(vals))


def group_delay_ms(frequencies_hz: Sequence[float], transfer: Sequence[complex]) -> np.ndarray:
    """Numerical group delay -d phase / d omega, in milliseconds."""
    f = np.asarray(frequencies_hz, dtype=float)
    h = np.asarray(transfer, dtype=np.complex128)
    if len(f) < 3:
        raise ValueError("need at least three frequencies for group delay")
    phase = np.unwrap(np.angle(h))
    omega = 2.0 * np.pi * f
    tau_s = -np.gradient(phase, omega)
    return 1e3 * tau_s
