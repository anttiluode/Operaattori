"""Gate-5 moving-operator mathematics.

A minimal matrix-valued dynamical system:

    h[t+1] = A[t] h[t] + (I - A[t]) 1 x[t]
    y_hat  = w[t]^T h[t]

with A[t] = diag(lambda[t]). The eigenvalues themselves are state:

    logit(lambda) = base + fast + slow

Prediction error is revealed only after prediction. A local sensitivity trace
p_i = d h_i / d logit(lambda_i) supplies an online eligibility signal; no BPTT
or replay buffer is used by the moving operator.

This module deliberately includes two worlds:
1. a smoothly drifting memory-kernel family, matched to this operator class;
2. an exact drifting-delay attack, where explicit context/attention should win.

The second world is load-bearing. Gate 5 is not allowed to become a
"moving matrices replace transformers" claim.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class MovingConfig:
    steps: int = 9000
    burnin: int = 800
    seeds: int = 8
    mode_taus: tuple[float, float] = (4.0, 16.0)
    fast_rate: float = 4.0
    slow_rate: float = 0.08
    fast_retention: float = 0.995
    readout_rate: float = 0.10
    context_kernel: int = 64
    context_delay: int = 32
    window_kernel_rate: float = 1.0
    window_delay_rate: float = 0.8
    attention_retention: float = 0.95
    attention_beta: float = 12.0

    def as_dict(self) -> dict:
        return asdict(self)


def _interp_anchors(steps: int, anchors: Iterable[float]) -> np.ndarray:
    anchors = np.asarray(tuple(anchors), dtype=float)
    return np.interp(
        np.arange(steps, dtype=float),
        np.linspace(0.0, steps - 1.0, len(anchors)),
        anchors,
    )


def drifting_kernel_world(
    seed: int,
    config: MovingConfig | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    c = config or MovingConfig()
    rng = np.random.default_rng(seed)
    x = rng.normal(size=c.steps)
    tau_fast = _interp_anchors(
        c.steps, (2.5, 4.0, 7.0, 12.0, 7.0, 4.0, 2.5)
    )
    tau_slow = _interp_anchors(
        c.steps, (9.0, 14.0, 24.0, 40.0, 24.0, 14.0, 9.0)
    )
    teacher_effective_tau = 0.65 * tau_fast + 0.35 * tau_slow
    state = np.zeros(2, dtype=float)
    y = np.zeros(c.steps, dtype=float)
    for t in range(c.steps):
        lam = np.exp(-1.0 / np.asarray([tau_fast[t], tau_slow[t]]))
        state = lam * state + (1.0 - lam) * x[t]
        y[t] = 0.65 * state[0] + 0.35 * state[1]
    return x, y, teacher_effective_tau


def drifting_delay_world(
    seed: int,
    config: MovingConfig | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    c = config or MovingConfig()
    max_lag = c.context_delay
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=c.steps + max_lag + 2)
    lags = np.rint(
        _interp_anchors(
            c.steps, (4.0, 6.0, 9.0, 14.0, 20.0, 14.0, 9.0, 6.0, 4.0)
        )
    ).astype(int)
    x = raw[max_lag : max_lag + c.steps]
    y = np.asarray(
        [raw[max_lag + t - int(lags[t])] for t in range(c.steps)],
        dtype=float,
    )
    return x, y, lags.astype(float)


class MovingDiagonalOperator:
    """Two-mode moving recurrent matrix with fast/slow operator state."""

    def __init__(
        self,
        config: MovingConfig,
        *,
        move_operator: bool = True,
    ) -> None:
        self.config = config
        taus = np.asarray(config.mode_taus, dtype=float)
        lam = np.exp(-1.0 / taus)
        self.base = np.log(lam / (1.0 - lam))
        self.fast = np.zeros_like(taus)
        self.slow = np.zeros_like(taus)
        self.h = np.zeros_like(taus)
        self.eligibility = np.zeros_like(taus)
        self.w = np.ones_like(taus) / len(taus)
        self.move_operator = move_operator

    def lambdas(self) -> np.ndarray:
        logits = np.clip(self.base + self.fast + self.slow, -3.5, 5.0)
        return 1.0 / (1.0 + np.exp(-logits))

    def operator_matrix(self) -> np.ndarray:
        return np.diag(self.lambdas())

    def mode_taus(self) -> np.ndarray:
        lam = np.clip(self.lambdas(), 1e-8, 0.99999)
        return -1.0 / np.log(lam)

    def effective_tau(self) -> float:
        weights = np.abs(self.w)
        return float(
            np.sum(weights * self.mode_taus())
            / (float(np.sum(weights)) + 1e-12)
        )

    def predict(self, x: float) -> float:
        lam = self.lambdas()
        h_prev = self.h.copy()
        dlam = lam * (1.0 - lam)
        self.eligibility = (
            lam * self.eligibility + dlam * (h_prev - x)
        )
        self.h = lam * h_prev + (1.0 - lam) * x
        return float(np.clip(self.w @ self.h, -5.0, 5.0))

    def learn(self, target: float, prediction: float) -> None:
        c = self.config
        error = float(np.clip(target - prediction, -3.0, 3.0))
        self.w += (
            c.readout_rate
            * error
            * self.h
            / (0.05 + float(self.h @ self.h))
        )
        norm = float(np.linalg.norm(self.w))
        if norm > 4.0:
            self.w *= 4.0 / norm
        if not self.move_operator:
            return
        local_credit = np.clip(
            error * self.w * self.eligibility,
            -1.0,
            1.0,
        )
        self.fast = (
            c.fast_retention * self.fast
            + c.fast_rate * local_credit
        )
        self.slow += c.slow_rate * local_credit
        total = np.clip(self.base + self.fast + self.slow, -3.5, 5.0)
        self.fast += total - (self.base + self.fast + self.slow)


class ExplicitContextLMS:
    """Strong scalar explicit-context / linear-attention ruler."""

    def __init__(self, window: int, rate: float) -> None:
        self.window = int(window)
        self.rate = float(rate)
        self.buffer: list[float] = []
        self.w = np.zeros(self.window, dtype=float)

    def _context(self) -> np.ndarray:
        v = np.zeros(self.window, dtype=float)
        n = min(self.window, len(self.buffer))
        if n:
            v[:n] = np.asarray(self.buffer[-n:][::-1], dtype=float)
        return v

    def predict(self, x: float) -> tuple[float, np.ndarray]:
        self.buffer.append(float(x))
        v = self._context()
        return float(self.w @ v), v

    def learn(
        self,
        target: float,
        prediction: float,
        context: np.ndarray,
    ) -> None:
        error = target - prediction
        self.w += (
            self.rate
            * error
            * context
            / (1e-6 + float(context @ context))
        )


class CausalLagAttention:
    """Transformer-like explicit-token attention over lag positions."""

    def __init__(
        self,
        window: int,
        retention: float,
        beta: float,
    ) -> None:
        self.window = int(window)
        self.retention = float(retention)
        self.beta = float(beta)
        self.buffer: list[float] = []
        self.score = np.zeros(self.window, dtype=float)

    def weights(self) -> np.ndarray:
        z = self.beta * (self.score - float(np.max(self.score)))
        a = np.exp(z)
        return a / (float(np.sum(a)) + 1e-12)

    def _lags(self) -> np.ndarray:
        v = np.zeros(self.window, dtype=float)
        available = min(self.window, max(0, len(self.buffer) - 1))
        for lag in range(1, available + 1):
            v[lag - 1] = self.buffer[-1 - lag]
        return v

    def predict(self, x: float) -> tuple[float, np.ndarray]:
        self.buffer.append(float(x))
        v = self._lags()
        return float(self.weights() @ v), v

    def learn(
        self,
        target: float,
        prediction: float,
        context: np.ndarray,
    ) -> None:
        del prediction
        self.score = (
            self.retention * self.score
            + (1.0 - self.retention) * context * target
        )

    def estimated_lag(self) -> float:
        lags = np.arange(1, self.window + 1, dtype=float)
        return float(self.weights() @ lags)


def _mse(pred: np.ndarray, target: np.ndarray, burnin: int) -> float:
    return float(np.mean((pred[burnin:] - target[burnin:]) ** 2))


def run_kernel_seed(
    seed: int,
    config: MovingConfig | None = None,
    *,
    keep_trace: bool = False,
) -> dict:
    c = config or MovingConfig()
    x, y, teacher_tau = drifting_kernel_world(seed, c)
    moving = MovingDiagonalOperator(c, move_operator=True)
    frozen = MovingDiagonalOperator(c, move_operator=False)
    context = ExplicitContextLMS(
        c.context_kernel, c.window_kernel_rate
    )
    pred_m = np.zeros(c.steps)
    pred_f = np.zeros(c.steps)
    pred_c = np.zeros(c.steps)
    learned_tau = np.zeros(c.steps)
    for t in range(c.steps):
        pm = moving.predict(float(x[t]))
        pf = frozen.predict(float(x[t]))
        pc, vc = context.predict(float(x[t]))
        pred_m[t] = pm
        pred_f[t] = pf
        pred_c[t] = pc
        moving.learn(float(y[t]), pm)
        frozen.learn(float(y[t]), pf)
        context.learn(float(y[t]), pc, vc)
        learned_tau[t] = moving.effective_tau()
    active = slice(c.burnin, None)
    moving_mse = _mse(pred_m, y, c.burnin)
    frozen_mse = _mse(pred_f, y, c.burnin)
    context_mse = _mse(pred_c, y, c.burnin)
    result = {
        "moving_mse": moving_mse,
        "frozen_mse": frozen_mse,
        "explicit_context_mse": context_mse,
        "moving_vs_frozen_ratio":
            moving_mse / (frozen_mse + 1e-12),
        "moving_vs_context_ratio":
            moving_mse / (context_mse + 1e-12),
        "operator_timescale_correlation": float(
            np.corrcoef(
                learned_tau[active],
                teacher_tau[active],
            )[0, 1]
        ),
        "moving_mutable_scalars": 10,
        "explicit_context_stored_scalars": 2 * c.context_kernel,
    }
    if keep_trace:
        stride = 50
        result["trace"] = {
            "step": list(range(0, c.steps, stride)),
            "teacher_effective_tau":
                teacher_tau[::stride].tolist(),
            "learned_effective_tau":
                learned_tau[::stride].tolist(),
        }
    return result


def run_delay_seed(
    seed: int,
    config: MovingConfig | None = None,
    *,
    keep_trace: bool = False,
) -> dict:
    c = config or MovingConfig()
    x, y, true_lag = drifting_delay_world(seed, c)
    moving = MovingDiagonalOperator(c, move_operator=True)
    frozen = MovingDiagonalOperator(c, move_operator=False)
    context = ExplicitContextLMS(
        c.context_delay, c.window_delay_rate
    )
    attention = CausalLagAttention(
        c.context_delay,
        c.attention_retention,
        c.attention_beta,
    )
    pred_m = np.zeros(c.steps)
    pred_f = np.zeros(c.steps)
    pred_c = np.zeros(c.steps)
    pred_a = np.zeros(c.steps)
    learned_tau = np.zeros(c.steps)
    attention_lag = np.zeros(c.steps)
    for t in range(c.steps):
        pm = moving.predict(float(x[t]))
        pf = frozen.predict(float(x[t]))
        pc, vc = context.predict(float(x[t]))
        pa, va = attention.predict(float(x[t]))
        pred_m[t] = pm
        pred_f[t] = pf
        pred_c[t] = pc
        pred_a[t] = pa
        moving.learn(float(y[t]), pm)
        frozen.learn(float(y[t]), pf)
        context.learn(float(y[t]), pc, vc)
        attention.learn(float(y[t]), pa, va)
        learned_tau[t] = moving.effective_tau()
        attention_lag[t] = attention.estimated_lag()
    active = slice(c.burnin, None)
    result = {
        "moving_mse": _mse(pred_m, y, c.burnin),
        "frozen_mse": _mse(pred_f, y, c.burnin),
        "explicit_context_mse": _mse(pred_c, y, c.burnin),
        "attention_mse": _mse(pred_a, y, c.burnin),
        "moving_lag_correlation": float(
            np.corrcoef(
                learned_tau[active],
                true_lag[active],
            )[0, 1]
        ),
        "attention_lag_correlation": float(
            np.corrcoef(
                attention_lag[active],
                true_lag[active],
            )[0, 1]
        ),
        "moving_mutable_scalars": 10,
        "attention_stored_scalars": 2 * c.context_delay,
    }
    if keep_trace:
        stride = 50
        result["trace"] = {
            "step": list(range(0, c.steps, stride)),
            "true_lag": true_lag[::stride].tolist(),
            "moving_effective_tau":
                learned_tau[::stride].tolist(),
            "attention_estimated_lag":
                attention_lag[::stride].tolist(),
        }
    return result


def _aggregate(rows: list[dict], keys: Iterable[str]) -> dict:
    out = {}
    for key in keys:
        values = np.asarray([float(row[key]) for row in rows])
        out[key] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
        }
    return out


def run_battery(config: MovingConfig | None = None) -> dict:
    c = config or MovingConfig()
    kernel_rows = [
        run_kernel_seed(seed, c, keep_trace=(seed == 0))
        for seed in range(c.seeds)
    ]
    delay_rows = [
        run_delay_seed(seed, c, keep_trace=(seed == 0))
        for seed in range(c.seeds)
    ]
    return {
        "config": c.as_dict(),
        "kernel": {
            "summary": _aggregate(
                kernel_rows,
                (
                    "moving_mse",
                    "frozen_mse",
                    "explicit_context_mse",
                    "moving_vs_frozen_ratio",
                    "moving_vs_context_ratio",
                    "operator_timescale_correlation",
                ),
            ),
            "seed0_trace": kernel_rows[0]["trace"],
            "moving_mutable_scalars": 10,
            "explicit_context_stored_scalars":
                2 * c.context_kernel,
        },
        "delay_attack": {
            "summary": _aggregate(
                delay_rows,
                (
                    "moving_mse",
                    "frozen_mse",
                    "explicit_context_mse",
                    "attention_mse",
                    "moving_lag_correlation",
                    "attention_lag_correlation",
                ),
            ),
            "seed0_trace": delay_rows[0]["trace"],
            "moving_mutable_scalars": 10,
            "attention_stored_scalars":
                2 * c.context_delay,
        },
    }


def write_json(path: str | Path, payload: dict) -> None:
    import json
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
