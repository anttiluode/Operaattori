"""Null models for auditing whether a real dendritic cable order is special.

Gate 11 established an order effect: reordering the same heterogeneous cable
segments changes transfer.  Gate 12 asks the harder question:

    is the biological ordering unusual relative to plausible reorderings,
    or is most of the effect explained by gross taper / endpoints?

Nothing here assigns a task objective to the neuron.  "Unusual" means only
that the real transfer signature lies in the tail of a permutation null.
It does not mean "better" or "optimized".
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cable_path import PassiveCableParams, cable_abcd


@dataclass(frozen=True)
class OrderAuditConfig:
    frequencies_hz: tuple[float, ...] = (1.0, 5.0, 15.0, 40.0, 100.0, 300.0)
    endpoint_fraction: float = 0.10
    coarse_windows_um: tuple[float, ...] = (10.0, 25.0, 50.0, 100.0)


def segment_midpoints_um(lengths_um: np.ndarray) -> np.ndarray:
    lengths = np.asarray(lengths_um, dtype=float)
    starts = np.concatenate([[0.0], np.cumsum(lengths[:-1])])
    return starts + 0.5 * lengths


def full_permutations(n: int, count: int, rng: np.random.Generator) -> np.ndarray:
    return np.stack([rng.permutation(n) for _ in range(count)], axis=0)


def endpoint_preserving_permutations(
    lengths_um: np.ndarray,
    count: int,
    rng: np.random.Generator,
    fraction: float = 0.10,
) -> np.ndarray:
    """Shuffle only the physical middle; keep proximal/distal endpoint zones."""
    lengths = np.asarray(lengths_um, dtype=float)
    n = len(lengths)
    mids = segment_midpoints_um(lengths)
    total = float(np.sum(lengths))
    fixed = (mids <= fraction * total) | (mids >= (1.0 - fraction) * total)
    middle = np.flatnonzero(~fixed)

    out = np.tile(np.arange(n, dtype=int), (count, 1))
    for row in out:
        row[middle] = rng.permutation(middle)
    return out


def within_window_permutations(
    lengths_um: np.ndarray,
    window_um: float,
    count: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Destroy fine order only within fixed proximal->distal distance windows.

    Every segment stays inside its original coarse spatial window, so the
    large-scale taper profile is preserved while fine sequence is randomized.
    """
    lengths = np.asarray(lengths_um, dtype=float)
    n = len(lengths)
    mids = segment_midpoints_um(lengths)
    bins = np.floor(mids / float(window_um)).astype(int)

    groups = [np.flatnonzero(bins == b) for b in np.unique(bins)]
    out = np.tile(np.arange(n, dtype=int), (count, 1))
    for row in out:
        for g in groups:
            if len(g) > 1:
                row[g] = rng.permutation(g)
    return out


def monotone_radius_orders(radii_um: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    radii = np.asarray(radii_um, dtype=float)
    # stable sort keeps original order among exact radius ties
    thick_to_thin = np.argsort(-radii, kind="stable")
    thin_to_thick = np.argsort(radii, kind="stable")
    return thick_to_thin.astype(int), thin_to_thick.astype(int)


def precompute_segment_matrices(
    lengths_um: np.ndarray,
    radii_um: np.ndarray,
    frequencies_hz: np.ndarray,
    params: PassiveCableParams | None = None,
) -> np.ndarray:
    lengths = np.asarray(lengths_um, dtype=float)
    radii = np.asarray(radii_um, dtype=float)
    freqs = np.asarray(frequencies_hz, dtype=float)
    if lengths.shape != radii.shape:
        raise ValueError("lengths and radii must have the same shape")

    out = np.empty((len(freqs), len(lengths), 2, 2), dtype=np.complex128)
    for fi, f in enumerate(freqs):
        for i, (l, r) in enumerate(zip(lengths, radii)):
            out[fi, i] = cable_abcd(float(l), float(r), float(f), params)
    return out


def compose_orders(segment_matrices: np.ndarray, orders: np.ndarray) -> np.ndarray:
    """Vectorized composition for many segment orders.

    Parameters
    ----------
    segment_matrices : F x N x 2 x 2
    orders           : P x N

    Returns
    -------
    P x F x 2 x 2
    """
    mats = np.asarray(segment_matrices, dtype=np.complex128)
    orders = np.asarray(orders, dtype=int)
    if mats.ndim != 4 or mats.shape[2:] != (2, 2):
        raise ValueError("segment_matrices must be F x N x 2 x 2")
    if orders.ndim == 1:
        orders = orders[None, :]
    if orders.ndim != 2 or orders.shape[1] != mats.shape[1]:
        raise ValueError("orders must be P x N")

    p = len(orders)
    f = mats.shape[0]
    out = np.broadcast_to(np.eye(2, dtype=np.complex128), (p, f, 2, 2)).copy()

    # Loop over path depth, vectorize over both permutations and frequencies.
    for j in range(orders.shape[1]):
        # selected: P x F x 2 x 2
        selected = np.transpose(mats[:, orders[:, j], :, :], (1, 0, 2, 3))
        out = out @ selected
    return out


def transfer_features(composed: np.ndarray) -> np.ndarray:
    """Return [gain dB, log10|Zin|, unwrapped gain phase] feature vectors."""
    M = np.asarray(composed, dtype=np.complex128)
    if M.ndim == 3:
        M = M[None, ...]
    if M.ndim != 4:
        raise ValueError("composed must be P x F x 2 x 2")

    A = M[:, :, 0, 0]
    C = M[:, :, 1, 0]
    gain = 1.0 / A
    zin = A / C

    gain_db = 20.0 * np.log10(np.maximum(np.abs(gain), 1e-300))
    log_z = np.log10(np.maximum(np.abs(zin), 1e-300))
    phase = np.unwrap(np.angle(gain), axis=1)
    return np.concatenate([gain_db, log_z, phase], axis=1)


def null_standardization(null_features: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = np.asarray(null_features, dtype=float)
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0)
    # Exclude dimensions with effectively no null variability rather than
    # manufacturing enormous z scores from floating-point crumbs.
    scale = np.maximum(np.median(std[std > 0]) if np.any(std > 0) else 0.0, 1.0)
    active = std > (1e-10 * scale)
    return mean, std, active


def standardized_rms(
    features: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    active: np.ndarray,
) -> np.ndarray:
    X = np.asarray(features, dtype=float)
    if X.ndim == 1:
        X = X[None, :]
    if not np.any(active):
        return np.zeros(len(X), dtype=float)
    z = (X[:, active] - mean[active]) / std[active]
    return np.sqrt(np.mean(z * z, axis=1))


def empirical_unusualness(
    real_features: np.ndarray,
    null_features: np.ndarray,
) -> dict[str, float]:
    """Tail probability for distance of real transfer signature from null mean."""
    null = np.asarray(null_features, dtype=float)
    real = np.asarray(real_features, dtype=float).reshape(1, -1)
    mean, std, active = null_standardization(null)

    null_scores = standardized_rms(null, mean, std, active)
    real_score = float(standardized_rms(real, mean, std, active)[0])

    # Upper-tail Monte-Carlo p: how many null orders are at least as unusual?
    tail_p = float((1 + np.sum(null_scores >= real_score)) / (len(null_scores) + 1))
    percentile = 1.0 - tail_p

    return {
        "real_standardized_rms": real_score,
        "null_score_mean": float(np.mean(null_scores)),
        "null_score_std": float(np.std(null_scores)),
        "empirical_upper_tail_p": tail_p,
        "unusualness_percentile": percentile,
        "active_feature_dimensions": int(np.sum(active)),
    }


def standardized_distance(
    a: np.ndarray,
    b: np.ndarray,
    null_features: np.ndarray,
) -> float:
    null = np.asarray(null_features, dtype=float)
    mean, std, active = null_standardization(null)
    if not np.any(active):
        return 0.0
    da = np.asarray(a, dtype=float)[active]
    db = np.asarray(b, dtype=float)[active]
    return float(np.sqrt(np.mean(((da - db) / std[active]) ** 2)))


def linear_taper_r2(lengths_um: np.ndarray, radii_um: np.ndarray) -> tuple[float, float]:
    """Simple gross-taper descriptor: radius ~ beta0 + beta1 * normalized distance."""
    mids = segment_midpoints_um(np.asarray(lengths_um, dtype=float))
    total = max(float(np.sum(lengths_um)), 1e-12)
    x = mids / total
    y = np.asarray(radii_um, dtype=float)
    X = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return float(beta[1]), float(r2)
