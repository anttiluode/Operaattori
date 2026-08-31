"""Gate-6 out-of-family moving-operator audit.

Gate 5 used a hidden two-mode drifting teacher and a two-mode moving student.
That was useful but dangerously matched. Gate 6 deliberately breaks the match:
the teacher is a six-mode positive memory kernel whose time constants and mode
weights drift independently.

The question is no longer "can the student rediscover its own family?" It is:

    can a lower-dimensional moving operator remain on the online-state
    budget/error frontier when the hidden world has richer temporal structure?

Explicit-context LMS rulers are swept across context budgets. Frozen operators
with the same recurrent capacity are included so any gain from motion is
visible separately from mere extra state.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class RichConfig:
    steps: int = 9000
    burnin: int = 900
    seeds: int = 6

    teacher_taus: tuple[float, ...] = (2.5, 5.0, 10.0, 20.0, 40.0, 80.0)
    moving_sizes: tuple[int, ...] = (2, 4, 8)
    context_windows: tuple[int, ...] = (4, 8, 16, 32, 64, 128)

    fast_rate: float = 3.0
    slow_rate: float = 0.055
    fast_retention: float = 0.996
    readout_rate: float = 0.075
    context_rate: float = 0.85

    def as_dict(self) -> dict:
        return asdict(self)


def _interp(steps: int, anchors: Iterable[float]) -> np.ndarray:
    anchors = np.asarray(tuple(anchors), dtype=float)
    return np.interp(
        np.arange(steps, dtype=float),
        np.linspace(0.0, steps - 1.0, len(anchors)),
        anchors,
    )


def rich_drifting_world(
    seed: int,
    config: RichConfig | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Six-mode teacher with independently drifting modes and mixture weights."""
    c = config or RichConfig()
    rng = np.random.default_rng(seed)
    x = rng.normal(size=c.steps)

    base = np.asarray(c.teacher_taus, dtype=float)
    n = len(base)
    t = np.linspace(0.0, 1.0, c.steps)

    global_scale = _interp(
        c.steps, (0.70, 0.95, 1.45, 1.75, 1.05, 0.78, 1.30, 0.72)
    )

    # Each mode drifts around the shared developmental trend with a different
    # phase. This breaks the "all eigenvalues scale together" shortcut.
    tau = np.empty((c.steps, n), dtype=float)
    for i in range(n):
        phase = 0.37 * i
        local = np.exp(
            0.18 * np.sin(2.0 * np.pi * (1.3 * t + phase))
            + 0.08 * np.sin(2.0 * np.pi * (3.1 * t + 0.19 * i))
        )
        tau[:, i] = base[i] * global_scale * local

    # Positive mixture weights wander among modes. Positive weights keep the
    # "effective memory horizon" interpretable as a weighted timescale.
    logits = np.empty((c.steps, n), dtype=float)
    for i in range(n):
        logits[:, i] = (
            0.55 * np.sin(2.0 * np.pi * (0.75 * t + 0.23 * i))
            + 0.35 * np.cos(2.0 * np.pi * (1.65 * t + 0.11 * i))
        )
    logits -= np.max(logits, axis=1, keepdims=True)
    weights = np.exp(logits)
    weights /= np.sum(weights, axis=1, keepdims=True)

    state = np.zeros(n, dtype=float)
    y = np.zeros(c.steps, dtype=float)
    horizon = np.zeros(c.steps, dtype=float)

    for k in range(c.steps):
        lam = np.exp(-1.0 / tau[k])
        state = lam * state + (1.0 - lam) * x[k]
        y[k] = float(weights[k] @ state)
        horizon[k] = float(weights[k] @ tau[k])

    return x, y, horizon


class BudgetMovingOperator:
    """Variable-width diagonal moving operator with local eligibility."""

    def __init__(
        self,
        modes: int,
        config: RichConfig,
        *,
        move_operator: bool,
    ) -> None:
        self.modes = int(modes)
        self.config = config
        initial_tau = np.geomspace(2.5, 80.0, self.modes)
        lam = np.exp(-1.0 / initial_tau)
        self.base = np.log(lam / (1.0 - lam))

        self.fast = np.zeros(self.modes, dtype=float)
        self.slow = np.zeros(self.modes, dtype=float)
        self.h = np.zeros(self.modes, dtype=float)
        self.eligibility = np.zeros(self.modes, dtype=float)
        self.w = np.ones(self.modes, dtype=float) / self.modes
        self.move_operator = bool(move_operator)

    @property
    def online_scalars(self) -> int:
        # h, eligibility, fast, slow, readout
        return 5 * self.modes

    def lambdas(self) -> np.ndarray:
        z = np.clip(self.base + self.fast + self.slow, -3.5, 5.2)
        return 1.0 / (1.0 + np.exp(-z))

    def mode_taus(self) -> np.ndarray:
        lam = np.clip(self.lambdas(), 1e-8, 0.99999)
        return -1.0 / np.log(lam)

    def effective_tau(self) -> float:
        a = np.abs(self.w)
        return float(
            a @ self.mode_taus() / (float(np.sum(a)) + 1e-12)
        )

    def predict(self, x: float) -> float:
        lam = self.lambdas()
        h0 = self.h.copy()
        dlam = lam * (1.0 - lam)

        self.eligibility = (
            lam * self.eligibility + dlam * (h0 - x)
        )
        self.h = lam * h0 + (1.0 - lam) * x
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
        # Scale structural step mildly with width so the wider operator does
        # not receive a hidden larger total update merely for having more modes.
        scale = 2.0 / np.sqrt(float(self.modes))
        self.fast = (
            c.fast_retention * self.fast
            + scale * c.fast_rate * local_credit
        )
        self.slow += scale * c.slow_rate * local_credit

        total = np.clip(self.base + self.fast + self.slow, -3.5, 5.2)
        self.fast += total - (self.base + self.fast + self.slow)


class ContextLMS:
    """Explicit causal context with one adaptive coefficient per stored token."""

    def __init__(self, window: int, rate: float) -> None:
        self.window = int(window)
        self.rate = float(rate)
        self.buffer: list[float] = []
        self.w = np.zeros(self.window, dtype=float)

    @property
    def online_scalars(self) -> int:
        # Stored token values plus adaptive context weights.
        return 2 * self.window

    def predict(self, x: float) -> tuple[float, np.ndarray]:
        self.buffer.append(float(x))
        v = np.zeros(self.window, dtype=float)
        n = min(self.window, len(self.buffer))
        if n:
            v[:n] = np.asarray(self.buffer[-n:][::-1], dtype=float)
        return float(self.w @ v), v

    def learn(
        self,
        target: float,
        prediction: float,
        context: np.ndarray,
    ) -> None:
        error = float(target - prediction)
        self.w += (
            self.rate
            * error
            * context
            / (1e-6 + float(context @ context))
        )


def _mse(pred: np.ndarray, target: np.ndarray, burnin: int) -> float:
    return float(np.mean((pred[burnin:] - target[burnin:]) ** 2))


def run_seed(
    seed: int,
    config: RichConfig | None = None,
    *,
    keep_trace: bool = False,
) -> dict:
    c = config or RichConfig()
    x, y, teacher_horizon = rich_drifting_world(seed, c)

    moving = {
        n: BudgetMovingOperator(n, c, move_operator=True)
        for n in c.moving_sizes
    }
    frozen = {
        n: BudgetMovingOperator(n, c, move_operator=False)
        for n in c.moving_sizes
    }
    contexts = {
        w: ContextLMS(w, c.context_rate)
        for w in c.context_windows
    }

    pred_m = {n: np.zeros(c.steps) for n in c.moving_sizes}
    pred_f = {n: np.zeros(c.steps) for n in c.moving_sizes}
    pred_c = {w: np.zeros(c.steps) for w in c.context_windows}
    tau_m = {n: np.zeros(c.steps) for n in c.moving_sizes}

    for t in range(c.steps):
        for n, model in moving.items():
            p = model.predict(float(x[t]))
            pred_m[n][t] = p
            model.learn(float(y[t]), p)
            tau_m[n][t] = model.effective_tau()

        for n, model in frozen.items():
            p = model.predict(float(x[t]))
            pred_f[n][t] = p
            model.learn(float(y[t]), p)

        for w, model in contexts.items():
            p, v = model.predict(float(x[t]))
            pred_c[w][t] = p
            model.learn(float(y[t]), p, v)

    active = slice(c.burnin, None)

    result: dict[str, object] = {
        "moving": {},
        "frozen": {},
        "context": {},
    }

    for n in c.moving_sizes:
        result["moving"][str(n)] = {
            "mse": _mse(pred_m[n], y, c.burnin),
            "online_scalars": moving[n].online_scalars,
            "horizon_correlation": float(
                np.corrcoef(
                    tau_m[n][active],
                    teacher_horizon[active],
                )[0, 1]
            ),
        }
        result["frozen"][str(n)] = {
            "mse": _mse(pred_f[n], y, c.burnin),
            "online_scalars": frozen[n].online_scalars,
        }

    for w in c.context_windows:
        result["context"][str(w)] = {
            "mse": _mse(pred_c[w], y, c.burnin),
            "online_scalars": contexts[w].online_scalars,
        }

    if keep_trace:
        stride = 50
        result["trace"] = {
            "step": list(range(0, c.steps, stride)),
            "teacher_horizon": teacher_horizon[::stride].tolist(),
            "moving_horizon": {
                str(n): tau_m[n][::stride].tolist()
                for n in c.moving_sizes
            },
        }

    return result


def _agg(values: Iterable[float]) -> dict[str, float]:
    a = np.asarray(tuple(values), dtype=float)
    return {
        "mean": float(np.mean(a)),
        "std": float(np.std(a)),
    }


def _pareto(points: list[dict]) -> list[dict]:
    """Budget/error Pareto frontier, lower is better in both dimensions."""
    out = []
    for p in points:
        dominated = False
        for q in points:
            if q is p:
                continue
            if (
                q["online_scalars"] <= p["online_scalars"]
                and q["mse"] <= p["mse"]
                and (
                    q["online_scalars"] < p["online_scalars"]
                    or q["mse"] < p["mse"]
                )
            ):
                dominated = True
                break
        if not dominated:
            out.append(dict(p))
    return sorted(out, key=lambda row: row["online_scalars"])


def run_battery(config: RichConfig | None = None) -> dict:
    c = config or RichConfig()
    rows = [
        run_seed(seed, c, keep_trace=(seed == 0))
        for seed in range(c.seeds)
    ]

    summary: dict[str, object] = {
        "moving": {},
        "frozen": {},
        "context": {},
    }

    for n in c.moving_sizes:
        key = str(n)
        summary["moving"][key] = {
            "mse": _agg(row["moving"][key]["mse"] for row in rows),
            "horizon_correlation": _agg(
                row["moving"][key]["horizon_correlation"] for row in rows
            ),
            "online_scalars": rows[0]["moving"][key]["online_scalars"],
        }
        summary["frozen"][key] = {
            "mse": _agg(row["frozen"][key]["mse"] for row in rows),
            "online_scalars": rows[0]["frozen"][key]["online_scalars"],
        }

    for w in c.context_windows:
        key = str(w)
        summary["context"][key] = {
            "mse": _agg(row["context"][key]["mse"] for row in rows),
            "online_scalars": rows[0]["context"][key]["online_scalars"],
        }

    frontier_points = []
    for n in c.moving_sizes:
        row = summary["moving"][str(n)]
        frontier_points.append(
            {
                "name": f"moving-{n}",
                "online_scalars": row["online_scalars"],
                "mse": row["mse"]["mean"],
            }
        )
    for w in c.context_windows:
        row = summary["context"][str(w)]
        frontier_points.append(
            {
                "name": f"context-{w}",
                "online_scalars": row["online_scalars"],
                "mse": row["mse"]["mean"],
            }
        )

    return {
        "config": c.as_dict(),
        "summary": summary,
        "pareto_frontier": _pareto(frontier_points),
        "seed0_trace": rows[0]["trace"],
    }


def write_json(path: str | Path, payload: dict) -> None:
    import json

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
