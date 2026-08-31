"""Gate-7: move the basis, not only the eigenvalues.

Gates 5-6 moved diagonal recurrence eigenvalues in a fixed basis. That is only
half a moving matrix. Gate 7 uses a hidden two-dimensional operator

    A(t) = R(phi(t)) diag(lambda_1(t), lambda_2(t)) R(phi(t))^T

whose eigenvalues AND eigenvectors drift.

A full moving student has three dynamical operator coordinates:
two eigenvalue logits and one basis angle. Online forward sensitivities provide
local eligibility for those coordinates. A diagonal moving student receives
the same signals and consequence but is forbidden to rotate its basis.

The target is the hidden teacher's two-dimensional state itself. This is
system identification, not a language benchmark; explicit-context LMS and RLS
attackers are included as non-morphological alternatives.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class RotationConfig:
    steps: int = 9000
    burnin: int = 900
    seeds: int = 6

    fast_retention: float = 0.996
    fast_rate: float = 0.12
    slow_rate: float = 0.0025

    context_windows: tuple[int, ...] = (4, 8, 16, 32)
    context_rate: float = 0.95

    def as_dict(self) -> dict:
        return asdict(self)


def _interp(steps: int, anchors: Iterable[float]) -> np.ndarray:
    anchors = np.asarray(tuple(anchors), dtype=float)
    return np.interp(
        np.arange(steps, dtype=float),
        np.linspace(0.0, steps - 1.0, len(anchors)),
        anchors,
    )


def _rotation(phi: float) -> tuple[np.ndarray, np.ndarray]:
    c = float(np.cos(phi))
    s = float(np.sin(phi))
    R = np.asarray([[c, -s], [s, c]], dtype=float)
    dR = np.asarray([[-s, -c], [c, -s]], dtype=float)
    return R, dR


def operator_from_coordinates(
    logit_lambda: np.ndarray,
    phi: float,
) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    lam = 1.0 / (1.0 + np.exp(-np.asarray(logit_lambda, dtype=float)))
    R, dR = _rotation(phi)
    D = np.diag(lam)
    A = R @ D @ R.T

    dlam = lam * (1.0 - lam)
    dD0 = np.diag([dlam[0], 0.0])
    dD1 = np.diag([0.0, dlam[1]])
    dA0 = R @ dD0 @ R.T
    dA1 = R @ dD1 @ R.T
    dAphi = dR @ D @ R.T + R @ D @ dR.T
    return A, (dA0, dA1, dAphi)


def rotating_teacher(
    seed: int,
    config: RotationConfig | None = None,
) -> dict[str, np.ndarray]:
    c = config or RotationConfig()
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(c.steps, 2))

    tau1 = _interp(
        c.steps, (3.0, 4.5, 7.0, 10.0, 6.5, 4.0, 3.0)
    )
    tau2 = _interp(
        c.steps, (18.0, 26.0, 42.0, 58.0, 36.0, 23.0, 18.0)
    )
    phi = _interp(
        c.steps, (-0.90, -0.35, 0.55, 1.00, 0.20, -0.75, -0.90)
    )

    y = np.zeros((c.steps, 2), dtype=float)
    A_trace = np.zeros((c.steps, 2, 2), dtype=float)
    h = np.zeros(2, dtype=float)
    I = np.eye(2)

    for t in range(c.steps):
        lam = np.exp(-1.0 / np.asarray([tau1[t], tau2[t]]))
        logits = np.log(lam / (1.0 - lam))
        A, _ = operator_from_coordinates(logits, float(phi[t]))
        h = A @ h + (I - A) @ x[t]
        y[t] = h
        A_trace[t] = A

    return {
        "x": x,
        "y": y,
        "phi": phi,
        "tau1": tau1,
        "tau2": tau2,
        "A": A_trace,
    }


class FullMovingOperator:
    """Two-state symmetric operator with moving spectrum and basis."""

    def __init__(
        self,
        config: RotationConfig,
        *,
        move_basis: bool = True,
        move_operator: bool = True,
    ) -> None:
        self.config = config
        initial_tau = np.asarray([4.0, 28.0])
        lam = np.exp(-1.0 / initial_tau)
        base_logits = np.log(lam / (1.0 - lam))

        self.base = np.asarray([base_logits[0], base_logits[1], 0.0])
        self.fast = np.zeros(3, dtype=float)
        self.slow = np.zeros(3, dtype=float)
        self.h = np.zeros(2, dtype=float)

        self.move_basis = bool(move_basis)
        self.move_operator = bool(move_operator)
        self.n_params = 3 if self.move_basis else 2
        self.P = np.zeros((2, self.n_params), dtype=float)

    @property
    def online_scalars(self) -> int:
        # hidden state + forward sensitivities + fast/slow operator state.
        return 2 + 2 * self.n_params + 2 * self.n_params

    def coordinates(self) -> tuple[np.ndarray, float]:
        total = self.base + self.fast + self.slow
        logits = np.clip(total[:2], -3.5, 5.2)
        phi = float(np.clip(total[2], -1.45, 1.45))
        return logits, phi

    def matrix_and_derivatives(
        self,
    ) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
        logits, phi = self.coordinates()
        A, deriv = operator_from_coordinates(logits, phi)
        if self.move_basis:
            return A, deriv
        # Fixed-basis diagonal attacker. Its phi is held at zero and there is
        # no angle eligibility coordinate.
        A0, deriv0 = operator_from_coordinates(logits, 0.0)
        return A0, deriv0[:2]

    def basis_angle(self) -> float:
        if not self.move_basis:
            return 0.0
        return self.coordinates()[1]

    def predict(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        A, deriv = self.matrix_and_derivatives()
        h0 = self.h.copy()
        delta = h0 - x

        new_P = A @ self.P
        for j, dA in enumerate(deriv):
            new_P[:, j] += dA @ delta

        self.P = new_P
        self.h = x + A @ delta
        return self.h.copy()

    def learn(self, target: np.ndarray, prediction: np.ndarray) -> None:
        if not self.move_operator:
            return

        c = self.config
        error = np.asarray(target, dtype=float) - np.asarray(
            prediction, dtype=float
        )

        # Gate-7 v1 normalized each operator coordinate independently.
        # That failed the preregistered basis-alignment criterion because the
        # eigenvalue and angle sensitivities are not orthogonal. Use the tiny
        # coupled Gauss-Newton / sensitivity solve instead:
        #
        #   (P^T P + ridge I) delta = P^T error
        #
        # This remains online and local to the three operator coordinates; it
        # does not backpropagate through the full history.
        gram = self.P.T @ self.P + 0.02 * np.eye(self.n_params)
        credit = np.linalg.solve(gram, self.P.T @ error)
        credit = np.clip(credit, -1.0, 1.0)

        self.fast[: self.n_params] = (
            c.fast_retention * self.fast[: self.n_params]
            + c.fast_rate * credit
        )
        self.slow[: self.n_params] += c.slow_rate * credit

        # Project back into the safe operator coordinate region.
        total = self.base + self.fast + self.slow
        projected = total.copy()
        projected[:2] = np.clip(projected[:2], -3.5, 5.2)
        if self.move_basis:
            projected[2] = np.clip(projected[2], -1.45, 1.45)
        else:
            projected[2] = 0.0
            self.fast[2] = 0.0
            self.slow[2] = 0.0
        self.fast += projected - total

    def matrix(self) -> np.ndarray:
        return self.matrix_and_derivatives()[0]


class VectorContextLMS:
    """Explicit vector token history with adaptive linear readout."""

    def __init__(self, window: int, rate: float) -> None:
        self.window = int(window)
        self.rate = float(rate)
        self.buffer: list[np.ndarray] = []
        self.W = np.zeros((2, 2 * self.window), dtype=float)

    @property
    def online_scalars(self) -> int:
        # 2W stored token scalars + 4W adaptive weights.
        return 6 * self.window

    def predict(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        self.buffer.append(np.asarray(x, dtype=float).copy())
        v = np.zeros(2 * self.window, dtype=float)
        n = min(self.window, len(self.buffer))
        if n:
            recent = np.asarray(self.buffer[-n:][::-1], dtype=float).reshape(-1)
            v[: 2 * n] = recent
        return self.W @ v, v

    def learn(
        self,
        target: np.ndarray,
        prediction: np.ndarray,
        context: np.ndarray,
    ) -> None:
        error = np.asarray(target) - np.asarray(prediction)
        self.W += (
            self.rate
            * np.outer(error, context)
            / (1e-6 + float(context @ context))
        )


class VectorRLS:
    """Expensive full-covariance explicit-context attacker."""

    def __init__(
        self,
        window: int = 32,
        forgetting: float = 0.995,
        ridge: float = 0.25,
    ) -> None:
        self.window = int(window)
        self.dim = 2 * self.window
        self.forgetting = float(forgetting)
        self.buffer: list[np.ndarray] = []
        self.W = np.zeros((2, self.dim), dtype=float)
        self.P = np.eye(self.dim, dtype=float) / float(ridge)

    @property
    def online_scalars(self) -> int:
        # buffer + output weights + inverse covariance.
        return self.dim + 2 * self.dim + self.dim * self.dim

    def predict(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        self.buffer.append(np.asarray(x, dtype=float).copy())
        v = np.zeros(self.dim, dtype=float)
        n = min(self.window, len(self.buffer))
        if n:
            recent = np.asarray(self.buffer[-n:][::-1], dtype=float).reshape(-1)
            v[: 2 * n] = recent
        return self.W @ v, v

    def learn(
        self,
        target: np.ndarray,
        prediction: np.ndarray,
        context: np.ndarray,
    ) -> None:
        Pv = self.P @ context
        denom = self.forgetting + float(context @ Pv)
        gain = Pv / (denom + 1e-12)
        error = np.asarray(target) - np.asarray(prediction)
        self.W += np.outer(error, gain)
        self.P = (
            self.P - np.outer(gain, context) @ self.P
        ) / self.forgetting


def _mse(pred: np.ndarray, target: np.ndarray, burnin: int) -> float:
    e = pred[burnin:] - target[burnin:]
    return float(np.mean(e * e))


def _operator_error(
    estimated: np.ndarray,
    teacher: np.ndarray,
    burnin: int,
) -> float:
    num = np.linalg.norm(
        estimated[burnin:] - teacher[burnin:], axis=(1, 2)
    )
    den = np.linalg.norm(teacher[burnin:], axis=(1, 2)) + 1e-12
    return float(np.mean(num / den))


def run_seed(
    seed: int,
    config: RotationConfig | None = None,
    *,
    keep_trace: bool = False,
) -> dict:
    c = config or RotationConfig()
    world = rotating_teacher(seed, c)
    x = world["x"]
    y = world["y"]

    full = FullMovingOperator(c, move_basis=True, move_operator=True)
    diagonal = FullMovingOperator(c, move_basis=False, move_operator=True)
    frozen = FullMovingOperator(c, move_basis=True, move_operator=False)

    contexts = {
        w: VectorContextLMS(w, c.context_rate)
        for w in c.context_windows
    }
    rls = VectorRLS(window=32, forgetting=0.995, ridge=0.25)

    pred_full = np.zeros_like(y)
    pred_diag = np.zeros_like(y)
    pred_frozen = np.zeros_like(y)
    pred_context = {
        w: np.zeros_like(y) for w in c.context_windows
    }
    pred_rls = np.zeros_like(y)

    A_full = np.zeros_like(world["A"])
    A_diag = np.zeros_like(world["A"])
    phi_full = np.zeros(c.steps, dtype=float)

    for t in range(c.steps):
        pf = full.predict(x[t])
        pd = diagonal.predict(x[t])
        pz = frozen.predict(x[t])
        pred_full[t] = pf
        pred_diag[t] = pd
        pred_frozen[t] = pz

        full.learn(y[t], pf)
        diagonal.learn(y[t], pd)
        frozen.learn(y[t], pz)

        A_full[t] = full.matrix()
        A_diag[t] = diagonal.matrix()
        phi_full[t] = full.basis_angle()

        for w, model in contexts.items():
            pc, vc = model.predict(x[t])
            pred_context[w][t] = pc
            model.learn(y[t], pc, vc)

        pr, vr = rls.predict(x[t])
        pred_rls[t] = pr
        rls.learn(y[t], pr, vr)

    active = slice(c.burnin, None)
    angle_delta = phi_full[active] - world["phi"][active]
    basis_alignment = float(np.mean(np.cos(2.0 * angle_delta)))

    result: dict[str, object] = {
        "full": {
            "mse": _mse(pred_full, y, c.burnin),
            "operator_error": _operator_error(
                A_full, world["A"], c.burnin
            ),
            "basis_alignment": basis_alignment,
            "online_scalars": full.online_scalars,
        },
        "diagonal": {
            "mse": _mse(pred_diag, y, c.burnin),
            "operator_error": _operator_error(
                A_diag, world["A"], c.burnin
            ),
            "online_scalars": diagonal.online_scalars,
        },
        "frozen": {
            "mse": _mse(pred_frozen, y, c.burnin),
            "online_scalars": frozen.online_scalars,
        },
        "context": {},
        "rls32": {
            "mse": _mse(pred_rls, y, c.burnin),
            "online_scalars": rls.online_scalars,
        },
    }

    for w in c.context_windows:
        result["context"][str(w)] = {
            "mse": _mse(pred_context[w], y, c.burnin),
            "online_scalars": contexts[w].online_scalars,
        }

    if keep_trace:
        stride = 50
        result["trace"] = {
            "step": list(range(0, c.steps, stride)),
            "teacher_phi": world["phi"][::stride].tolist(),
            "learned_phi": phi_full[::stride].tolist(),
        }

    return result


def _agg(values: Iterable[float]) -> dict[str, float]:
    a = np.asarray(tuple(values), dtype=float)
    return {"mean": float(np.mean(a)), "std": float(np.std(a))}


def _pareto(points: list[dict]) -> list[dict]:
    out = []
    for p in points:
        if any(
            q is not p
            and q["online_scalars"] <= p["online_scalars"]
            and q["mse"] <= p["mse"]
            and (
                q["online_scalars"] < p["online_scalars"]
                or q["mse"] < p["mse"]
            )
            for q in points
        ):
            continue
        out.append(dict(p))
    return sorted(out, key=lambda row: row["online_scalars"])


def run_battery(config: RotationConfig | None = None) -> dict:
    c = config or RotationConfig()
    rows = [
        run_seed(seed, c, keep_trace=(seed == 0))
        for seed in range(c.seeds)
    ]

    summary = {
        "full": {
            "mse": _agg(row["full"]["mse"] for row in rows),
            "operator_error": _agg(
                row["full"]["operator_error"] for row in rows
            ),
            "basis_alignment": _agg(
                row["full"]["basis_alignment"] for row in rows
            ),
            "online_scalars": rows[0]["full"]["online_scalars"],
        },
        "diagonal": {
            "mse": _agg(row["diagonal"]["mse"] for row in rows),
            "operator_error": _agg(
                row["diagonal"]["operator_error"] for row in rows
            ),
            "online_scalars": rows[0]["diagonal"]["online_scalars"],
        },
        "frozen": {
            "mse": _agg(row["frozen"]["mse"] for row in rows),
            "online_scalars": rows[0]["frozen"]["online_scalars"],
        },
        "context": {},
        "rls32": {
            "mse": _agg(row["rls32"]["mse"] for row in rows),
            "online_scalars": rows[0]["rls32"]["online_scalars"],
        },
    }

    for w in c.context_windows:
        summary["context"][str(w)] = {
            "mse": _agg(row["context"][str(w)]["mse"] for row in rows),
            "online_scalars":
                rows[0]["context"][str(w)]["online_scalars"],
        }

    points = [
        {
            "name": "moving-full",
            "online_scalars": summary["full"]["online_scalars"],
            "mse": summary["full"]["mse"]["mean"],
        },
        {
            "name": "moving-diagonal",
            "online_scalars": summary["diagonal"]["online_scalars"],
            "mse": summary["diagonal"]["mse"]["mean"],
        },
    ]
    for w in c.context_windows:
        row = summary["context"][str(w)]
        points.append(
            {
                "name": f"context-{w}",
                "online_scalars": row["online_scalars"],
                "mse": row["mse"]["mean"],
            }
        )
    points.append(
        {
            "name": "rls-32",
            "online_scalars": summary["rls32"]["online_scalars"],
            "mse": summary["rls32"]["mse"]["mean"],
        }
    )

    return {
        "config": c.as_dict(),
        "summary": summary,
        "pareto_frontier": _pareto(points),
        "seed0_trace": rows[0]["trace"],
    }


def write_json(path: str | Path, payload: dict) -> None:
    import json

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
