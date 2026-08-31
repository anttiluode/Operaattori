"""Core mechanisms for the first Operaattori gates.

The first substrate is intentionally a one-dimensional material strip, not a
literal dendrite. It is small enough that every state variable and the exact
frozen linear operator can be inspected.

Development:
    input -> fast field -> local eligibility -> hysteretic material
                                      ^              |
                                      |              v
                                      +--- transport-+

After development, the fast and eligibility states are washed out. The
persistent material is then frozen and probed.

The primary histories are stronger than a simple AB-vs-BA test: they have the
same marginals, same first symbol and same last symbol.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

H0 = "AAABBB"
H1 = "ABBAAB"


@dataclass(frozen=True)
class Config:
    n: int = 32
    dt: float = 0.08

    prelude_steps: int = 80
    block_steps: int = 70
    gap_steps: int = 25
    common_suffix_steps: int = 120
    washout_steps: int = 800

    fast_decay: float = 0.22
    diffusion: float = 0.45
    input_gain: float = 1.8
    feedback: float = 1.4

    eligibility_decay: float = 0.12
    eligibility_gain: float = 0.40

    morph_gain: float = 0.90
    morph_decay: float = 0.03
    morph_threshold: float = 0.28
    bistable_gain: float = 0.65
    morph_noise: float = 0.003

    resource_recovery: float = 0.01
    resource_cost: float = 0.20

    # Frozen-operator probe. These are intentionally more conductive than the
    # developmental strip so that a one-cell structural difference has a
    # measurable downstream impulse response.
    probe_dt: float = 0.01
    probe_steps: int = 1200
    probe_pulse_steps: int = 60
    probe_diffusion: float = 0.80
    probe_decay: float = 0.12
    probe_base_conductance: float = 0.20
    probe_morph_conductance: float = 8.0
    probe_input_gain: float = 1.5

    @property
    def washout_eligibility_time_constants(self) -> float:
        return self.washout_steps * self.dt * self.eligibility_decay

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


def history_invariants(h0: str = H0, h1: str = H1) -> dict[str, bool | int]:
    return {
        "same_length": len(h0) == len(h1),
        "same_A_count": h0.count("A") == h1.count("A"),
        "same_B_count": h0.count("B") == h1.count("B"),
        "same_first": h0[0] == h1[0],
        "same_last": h0[-1] == h1[-1],
        "different_order": h0 != h1,
    }


def schedule(history: str, config: Config) -> np.ndarray:
    """Return the complete two-channel development protocol.

    Both histories receive the exact same common suffix and silence. The
    common suffix alternates A/B every step and ends on B because its length is
    even.
    """
    rows: list[tuple[float, float]] = [(0.0, 0.0)] * config.prelude_steps
    for symbol in history:
        if symbol == "A":
            stim = (1.0, 0.0)
        elif symbol == "B":
            stim = (0.0, 1.0)
        else:
            raise ValueError(f"unknown history symbol {symbol!r}")
        rows.extend([stim] * config.block_steps)
        rows.extend([(0.0, 0.0)] * config.gap_steps)

    rows.extend(
        [(1.0, 0.0) if i % 2 == 0 else (0.0, 1.0)
         for i in range(config.common_suffix_steps)]
    )
    rows.extend([(0.0, 0.0)] * config.washout_steps)
    return np.asarray(rows, dtype=float)


class MorphologicalSubstrate:
    """A local excitable strip with slow hysteretic material."""

    def __init__(
        self,
        config: Config,
        clone_seed: int,
        noise_seed: int,
    ) -> None:
        self.config = config
        init_rng = np.random.default_rng(clone_seed)
        self.noise_rng = np.random.default_rng(noise_seed)

        self.fast = np.zeros(config.n, dtype=float)
        self.eligibility = np.zeros(config.n, dtype=float)
        self.morphology = np.clip(
            0.04 + 0.005 * init_rng.normal(size=config.n),
            0.0,
            1.0,
        )
        self.resource = np.ones(config.n, dtype=float)

    def step(self, stimulus: Sequence[float]) -> None:
        c = self.config
        u = self.fast
        m = self.morphology

        # Local edge conductance. The material written by previous activity
        # already changes the next fast-field update.
        node_g = 0.18 + c.feedback * m
        edge_g = 0.5 * (node_g[:-1] + node_g[1:])
        flux = edge_g * (u[:-1] - u[1:])
        divergence = np.zeros_like(u)
        divergence[:-1] -= flux
        divergence[1:] += flux

        du = c.diffusion * divergence - c.fast_decay * u
        du[2] += c.input_gain * float(stimulus[0])
        du[-3] += c.input_gain * float(stimulus[1])
        self.fast = u + c.dt * du

        # Eligibility is driven by energy plus local spatial current. It is
        # deliberately labile and must be gone before the order readout.
        grad = np.zeros_like(u)
        grad[1:-1] = 0.5 * np.abs(
            self.fast[2:] - self.fast[:-2]
        )
        local_drive = np.tanh(
            2.0 * (self.fast * self.fast + 0.5 * grad)
        )
        self.eligibility += c.dt * (
            c.eligibility_gain * local_drive
            - c.eligibility_decay * self.eligibility
        )

        # Bistable material. Below threshold it tends to disappear; above it
        # the local high-material state can survive silence. Activity can push
        # a site across the threshold, and previously grown material changes
        # later activity through conductance above.
        bistable = (
            c.bistable_gain
            * m
            * (1.0 - m)
            * (m - c.morph_threshold)
        )
        growth = (
            c.morph_gain
            * self.eligibility
            * self.resource
            * (1.0 - m)
        )
        dm = bistable + growth - c.morph_decay * m
        dm += c.morph_noise * self.noise_rng.normal(size=c.n)
        new_m = np.clip(m + c.dt * dm, 0.0, 1.0)

        # Local finite resource. It is not used as memory: it recovers toward
        # one and only the persistent morphology is read at the order gate.
        consumed = np.maximum(new_m - m, 0.0)
        self.resource = np.clip(
            self.resource
            + c.dt * c.resource_recovery * (1.0 - self.resource)
            - c.resource_cost * consumed,
            0.0,
            1.0,
        )
        self.morphology = new_m

    def run(self, protocol: np.ndarray) -> None:
        for stimulus in protocol:
            self.step(stimulus)


def run_morphology(
    history: str,
    clone_seed: int,
    noise_seed: int,
    config: Config | None = None,
) -> MorphologicalSubstrate:
    c = config or Config()
    world = MorphologicalSubstrate(c, clone_seed, noise_seed)
    world.run(schedule(history, c))
    return world


class ContractiveMatrixRuler:
    """Same-dimensional stable linear state ruler.

    The matrix parameters are fixed across all clones. Clone/noise seeds change
    only initial microscopic state and process noise.
    """

    def __init__(
        self,
        config: Config,
        clone_seed: int,
        noise_seed: int,
        matrix_seed: int = 9137,
    ) -> None:
        self.config = config
        rng = np.random.default_rng(matrix_seed)
        q, _ = np.linalg.qr(rng.normal(size=(config.n, config.n)))
        retention = float(
            np.exp(-config.dt * config.eligibility_decay)
        )
        eig = retention * np.linspace(0.85, 1.0, config.n)
        self.A = q @ np.diag(eig) @ q.T
        self.B = rng.normal(size=(config.n, 2)) / np.sqrt(2.0)

        init_rng = np.random.default_rng(clone_seed)
        self.noise_rng = np.random.default_rng(noise_seed)
        self.state = 0.005 * init_rng.normal(size=config.n)

    def run(self, protocol: np.ndarray) -> np.ndarray:
        for stimulus in protocol:
            self.state = (
                self.A @ self.state
                + self.B @ stimulus
                + 0.003 * self.noise_rng.normal(size=self.config.n)
            )
        return self.state.copy()


class BistableStateRuler:
    """Same-capacity non-geometric hysteretic state attacker.

    This is supposed to be dangerous. If persistent hysteresis alone can store
    the order, this control should expose that.
    """

    def __init__(
        self,
        config: Config,
        clone_seed: int,
        noise_seed: int,
        ruler_seed: int = 811,
    ) -> None:
        self.config = config
        fixed_rng = np.random.default_rng(ruler_seed)
        self.B = fixed_rng.normal(size=(config.n, 2))

        init_rng = np.random.default_rng(clone_seed)
        self.noise_rng = np.random.default_rng(noise_seed)
        self.state = np.clip(
            0.04 + 0.005 * init_rng.normal(size=config.n),
            0.0,
            1.0,
        )

    def run(self, protocol: np.ndarray) -> np.ndarray:
        c = self.config
        for stimulus in protocol:
            drive = np.tanh(np.maximum(0.0, self.B @ stimulus))
            x = self.state
            bistable = (
                c.bistable_gain
                * x
                * (1.0 - x)
                * (x - c.morph_threshold)
            )
            growth = 0.25 * drive * (1.0 - x)
            dx = bistable + growth - c.morph_decay * x
            dx += c.morph_noise * self.noise_rng.normal(size=c.n)
            self.state = np.clip(x + c.dt * dx, 0.0, 1.0)
        return self.state.copy()


def nearest_centroid_loo_pairs(
    arm0: np.ndarray,
    arm1: np.ndarray,
) -> float:
    """Leave one clone pair out; Euclidean nearest-centroid classification."""
    if arm0.shape != arm1.shape:
        raise ValueError("matched arms must have the same shape")
    n = arm0.shape[0]
    correct = 0
    for i in range(n):
        keep = np.ones(n, dtype=bool)
        keep[i] = False
        c0 = np.mean(arm0[keep], axis=0)
        c1 = np.mean(arm1[keep], axis=0)
        d00 = float(np.sum((arm0[i] - c0) ** 2))
        d01 = float(np.sum((arm0[i] - c1) ** 2))
        d10 = float(np.sum((arm1[i] - c0) ** 2))
        d11 = float(np.sum((arm1[i] - c1) ** 2))
        correct += int(d00 <= d01)
        correct += int(d11 < d10)
    return correct / float(2 * n)


def unit_mass_shape(morphology: np.ndarray) -> np.ndarray:
    return morphology / (float(np.sum(morphology)) + 1e-12)


def mass_match_pair(
    a: np.ndarray,
    b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Scale a pair to identical total material without changing addresses."""
    target = 0.5 * (float(np.sum(a)) + float(np.sum(b)))
    aa = np.asarray(a, dtype=float) * target / (float(np.sum(a)) + 1e-12)
    bb = np.asarray(b, dtype=float) * target / (float(np.sum(b)) + 1e-12)
    return aa, bb


def _probe_laplacian(
    morphology: np.ndarray,
    config: Config,
) -> np.ndarray:
    n = config.n
    node_g = (
        config.probe_base_conductance
        + config.probe_morph_conductance * morphology
    )
    edge_g = 0.5 * (node_g[:-1] + node_g[1:])
    L = np.zeros((n, n), dtype=float)
    for i, ge in enumerate(edge_g):
        L[i, i] -= ge
        L[i, i + 1] += ge
        L[i + 1, i] += ge
        L[i + 1, i + 1] -= ge
    return L


def frozen_operator_matrix(
    morphology: np.ndarray,
    config: Config | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Exact frozen fast-state matrix A[m] and two-channel input matrix B."""
    c = config or Config()
    L = _probe_laplacian(np.asarray(morphology, dtype=float), c)
    A = np.eye(c.n) + c.probe_dt * (
        c.probe_diffusion * L - c.probe_decay * np.eye(c.n)
    )
    B = np.zeros((c.n, 2), dtype=float)
    B[2, 0] = c.probe_dt * c.probe_input_gain
    B[-3, 1] = c.probe_dt * c.probe_input_gain
    return A, B


def _probe_protocol(port: int, config: Config) -> np.ndarray:
    protocol = np.zeros((config.probe_steps, 2), dtype=float)
    protocol[: config.probe_pulse_steps, port] = 1.0
    return protocol


def probe_direct(
    morphology: np.ndarray,
    config: Config | None = None,
) -> np.ndarray:
    """Direct spatial simulation; returns A- and B-probe soma traces."""
    c = config or Config()
    L = _probe_laplacian(np.asarray(morphology, dtype=float), c)
    traces: list[np.ndarray] = []
    soma = c.n // 2
    for port in (0, 1):
        u = np.zeros(c.n, dtype=float)
        trace = np.zeros(c.probe_steps, dtype=float)
        for t, stimulus in enumerate(_probe_protocol(port, c)):
            du = c.probe_diffusion * (L @ u) - c.probe_decay * u
            du[2] += c.probe_input_gain * stimulus[0]
            du[-3] += c.probe_input_gain * stimulus[1]
            u = u + c.probe_dt * du
            trace[t] = u[soma]
        traces.append(trace)
    return np.concatenate(traces)


def probe_matrix(
    morphology: np.ndarray,
    config: Config | None = None,
) -> np.ndarray:
    """Exact matrix replay of the same frozen spatial probe."""
    c = config or Config()
    A, B = frozen_operator_matrix(morphology, c)
    traces: list[np.ndarray] = []
    soma = c.n // 2
    for port in (0, 1):
        u = np.zeros(c.n, dtype=float)
        trace = np.zeros(c.probe_steps, dtype=float)
        for t, stimulus in enumerate(_probe_protocol(port, c)):
            u = A @ u + B @ stimulus
            trace[t] = u[soma]
        traces.append(trace)
    return np.concatenate(traces)


def relative_l2(a: np.ndarray, b: np.ndarray) -> float:
    denom = 0.5 * (float(np.linalg.norm(a)) + float(np.linalg.norm(b)))
    return float(np.linalg.norm(a - b) / (denom + 1e-12))


def write_json(path: str | Path, payload: dict) -> None:
    import json

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
