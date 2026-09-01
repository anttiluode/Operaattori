"""Exact passive-cable transport on a branching tree.

Gate 11 treated root-to-tip paths as independent serial cable products. Gate 13
adds the missing neuron-specific operation: branch junctions.

A cable edge contributes a symmetric two-node admittance matrix. Descendant
subtrees are eliminated exactly by Schur complements on the tree. For a chosen
soma-to-tip route, every off-path subtree is equivalently a shunt admittance

    S(Y) = [[1, 0],
            [Y, 1]]

inserted between serial cable ABCD matrices.

This gives two independently computed views of the same passive tree:
1) whole-tree elimination;
2) cable-matrix / side-shunt composition along a selected route.

Their agreement is the representation receipt.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cable_path import PassiveCableParams, cable_abcd, cable_constants


@dataclass
class TreeFrequencyState:
    frequency_hz: float
    a_effective: np.ndarray
    transfer_to_clamp: np.ndarray
    y11: np.ndarray
    y12: np.ndarray
    y22: np.ndarray
    edge_active: np.ndarray


def cable_edge_admittance(
    length_um: float,
    radius_um: float,
    frequency_hz: float,
    params: PassiveCableParams | None = None,
) -> tuple[complex, complex, complex]:
    """Return (Y11, Y12, Y22) with currents entering the cable at both ends."""
    p = params or PassiveCableParams()
    length_cm = float(length_um) * 1e-4
    gamma, z0 = cable_constants(radius_um, frequency_hz, p)
    gl = gamma * length_cm
    sh = np.sinh(gl)
    ch = np.cosh(gl)
    if abs(sh) < 1e-18:
        raise ValueError("electrical segment too short for stable admittance form")
    y11 = ch / (z0 * sh)
    y12 = -1.0 / (z0 * sh)
    y22 = y11
    return complex(y11), complex(y12), complex(y22)


def active_child_counts(parents: np.ndarray, active: np.ndarray) -> np.ndarray:
    parents = np.asarray(parents, dtype=np.int64)
    active = np.asarray(active, dtype=bool)
    out = np.zeros(len(parents), dtype=np.int64)
    for i in range(1, len(parents)):
        if not active[i]:
            continue
        p = int(parents[i])
        if p >= 0 and active[p]:
            out[p] += 1
    return out


def solve_tree_frequency(
    parents: np.ndarray,
    lengths_um: np.ndarray,
    radii_um: np.ndarray,
    active: np.ndarray,
    clamped: np.ndarray,
    frequency_hz: float,
    params: PassiveCableParams | None = None,
) -> TreeFrequencyState:
    """Precompute exact current transfer from every active node to a clamped root.

    parents must be parent-before-child. Every active, non-clamped node must
    have an active parent and a physical cable edge stored at its own index.
    """
    parents = np.asarray(parents, dtype=np.int64)
    lengths = np.asarray(lengths_um, dtype=float)
    radii = np.asarray(radii_um, dtype=float)
    active = np.asarray(active, dtype=bool)
    clamped = np.asarray(clamped, dtype=bool)
    n = len(parents)

    y11 = np.zeros(n, dtype=np.complex128)
    y12 = np.zeros(n, dtype=np.complex128)
    y22 = np.zeros(n, dtype=np.complex128)
    edge_active = np.zeros(n, dtype=bool)

    for i in range(1, n):
        if not active[i] or clamped[i]:
            continue
        pidx = int(parents[i])
        if pidx < 0 or not active[pidx]:
            raise ValueError(f"active node {i} lacks active parent")
        vals = cable_edge_admittance(lengths[i], radii[i], frequency_hz, params)
        y11[i], y12[i], y22[i] = vals
        edge_active[i] = True

    # Base nodal diagonal after writing all incident cable Y11/Y22 terms.
    a = np.zeros(n, dtype=np.complex128)
    for i in range(1, n):
        if edge_active[i]:
            a[i] += y22[i]
            pidx = int(parents[i])
            if not clamped[pidx]:
                a[pidx] += y11[i]

    # Exact leaf-to-root Schur elimination of descendant voltages.
    for i in range(n - 1, 0, -1):
        if not edge_active[i]:
            continue
        if abs(a[i]) < 1e-300:
            raise FloatingPointError(f"zero effective admittance at node {i}")
        pidx = int(parents[i])
        if not clamped[pidx]:
            a[pidx] -= (y12[i] * y12[i]) / a[i]

    # Transfer from unit current injected at each node to its first clamped
    # ancestor. Sign follows the cable-port convention; magnitude/phase matter.
    T = np.zeros(n, dtype=np.complex128)
    for i in range(1, n):
        if not edge_active[i]:
            continue
        pidx = int(parents[i])
        if clamped[pidx]:
            T[i] = y12[i] / a[i]
        else:
            T[i] = T[pidx] * (-y12[i] / a[i])

    return TreeFrequencyState(
        frequency_hz=float(frequency_hz),
        a_effective=a,
        transfer_to_clamp=T,
        y11=y11,
        y12=y12,
        y22=y22,
        edge_active=edge_active,
    )


def path_to_clamp(
    parents: np.ndarray,
    clamped: np.ndarray,
    node: int,
) -> list[int]:
    """Return physical edge-child indices from proximal to distal."""
    parents = np.asarray(parents, dtype=np.int64)
    clamped = np.asarray(clamped, dtype=bool)
    path = []
    i = int(node)
    while not clamped[i]:
        path.append(i)
        i = int(parents[i])
        if i < 0:
            raise ValueError("node has no clamped ancestor")
    path.reverse()
    return path


def isolated_path_transfer(
    path_nodes: list[int],
    lengths_um: np.ndarray,
    radii_um: np.ndarray,
    frequency_hz: float,
    params: PassiveCableParams | None = None,
) -> complex:
    M = np.eye(2, dtype=np.complex128)
    for i in path_nodes:
        M = M @ cable_abcd(
            float(lengths_um[i]),
            float(radii_um[i]),
            float(frequency_hz),
            params,
        )
    return complex(-1.0 / M[0, 0])


def child_input_admittance_at_parent(
    child: int,
    state: TreeFrequencyState,
) -> complex:
    """Effective admittance of a child's complete subtree seen by its parent."""
    i = int(child)
    return complex(
        state.y11[i]
        - (state.y12[i] * state.y12[i]) / state.a_effective[i]
    )


def side_shunt_admittances(
    parents: np.ndarray,
    active: np.ndarray,
    clamped: np.ndarray,
    path_nodes: list[int],
    state: TreeFrequencyState,
) -> list[complex]:
    """Off-path subtree admittance attached at each node of a chosen path."""
    parents = np.asarray(parents, dtype=np.int64)
    active = np.asarray(active, dtype=bool)
    path_set_next = {
        path_nodes[j]: path_nodes[j + 1]
        for j in range(len(path_nodes) - 1)
    }
    children: dict[int, list[int]] = {}
    for i in range(1, len(parents)):
        if not state.edge_active[i]:
            continue
        children.setdefault(int(parents[i]), []).append(i)

    out = []
    for node in path_nodes:
        next_node = path_set_next.get(node)
        y = 0.0 + 0.0j
        for child in children.get(node, []):
            if child == next_node:
                continue
            y += child_input_admittance_at_parent(child, state)
        out.append(complex(y))
    return out


def shunt_abcd(admittance: complex) -> np.ndarray:
    return np.asarray(
        [[1.0 + 0.0j, 0.0 + 0.0j], [complex(admittance), 1.0 + 0.0j]],
        dtype=np.complex128,
    )


def path_with_side_shunts_transfer(
    path_nodes: list[int],
    lengths_um: np.ndarray,
    radii_um: np.ndarray,
    side_shunts: list[complex],
    frequency_hz: float,
    params: PassiveCableParams | None = None,
) -> complex:
    """Compose cable edges and exact off-path shunts along one route."""
    if len(path_nodes) != len(side_shunts):
        raise ValueError("one side-shunt value required per path node")

    M = np.eye(2, dtype=np.complex128)
    for j, node in enumerate(path_nodes):
        M = M @ cable_abcd(
            float(lengths_um[node]),
            float(radii_um[node]),
            float(frequency_hz),
            params,
        )
        # Side branches live at the distal node reached by this segment.
        if abs(side_shunts[j]) > 0:
            M = M @ shunt_abcd(side_shunts[j])
    return complex(-1.0 / M[0, 0])


def transfer_signature_features(transfer: np.ndarray) -> np.ndarray:
    """Per-node frequency signature: log magnitude + unwrapped phase."""
    H = np.asarray(transfer, dtype=np.complex128)
    mag_db = 20.0 * np.log10(np.maximum(np.abs(H), 1e-300))
    phase = np.unwrap(np.angle(H), axis=1)
    return np.concatenate([mag_db, phase], axis=1)


def standardized_effective_rank(features: np.ndarray) -> float:
    """Participation-ratio rank of centered, column-standardized signatures."""
    X = np.asarray(features, dtype=float)
    X = X - np.mean(X, axis=0, keepdims=True)
    std = np.std(X, axis=0)
    active = std > 1e-12
    if np.sum(active) == 0:
        return 0.0
    X = X[:, active] / std[active]
    s = np.linalg.svd(X, compute_uv=False)
    power = s * s
    denom = float(np.sum(power * power))
    if denom <= 0:
        return 0.0
    return float((np.sum(power) ** 2) / denom)
