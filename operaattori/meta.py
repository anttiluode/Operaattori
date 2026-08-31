"""Gate-8: learn the developmental law G_theta, not the final operator.

Outer loop:
    choose one small rule parameter vector theta on TRAIN worlds.

Inner loop on every world:
    start from the same operator;
    receive x_t;
    predict;
    reveal y_t;
    move A_t online using theta.

The selected theta is then frozen and evaluated on held-out TEST worlds.
No final operator state, readout, trajectory, or teacher parameters are carried
from train worlds into test worlds.

This is intentionally a tiny evolutionary/meta-learning analogue of
"genes are the rule that moves the operator."
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .rich import BudgetMovingOperator, ContextLMS, RichConfig


@dataclass(frozen=True)
class MetaConfig:
    steps: int = 4200
    burnin: int = 500
    modes: int = 4
    candidate_count: int = 40
    train_worlds: int = 8
    test_worlds: int = 8
    meta_seed: int = 20260831
    context_window: int = 64

    def as_dict(self) -> dict:
        return asdict(self)


def _interp(steps: int, anchors: np.ndarray) -> np.ndarray:
    return np.interp(
        np.arange(steps, dtype=float),
        np.linspace(0.0, steps - 1.0, len(anchors)),
        np.asarray(anchors, dtype=float),
    )


def meta_world(
    world_seed: int,
    config: MetaConfig | None = None,
) -> dict[str, np.ndarray]:
    """Random member of a six-mode smoothly drifting world family."""
    c = config or MetaConfig()
    rng = np.random.default_rng(world_seed)
    x = rng.normal(size=c.steps)

    base_tau = np.asarray([2.5, 5.5, 12.0, 26.0, 55.0, 105.0])
    n = len(base_tau)
    t = np.linspace(0.0, 1.0, c.steps)

    # Each world gets a different coarse developmental trajectory.
    anchors = np.exp(rng.uniform(np.log(0.62), np.log(1.70), size=7))
    global_scale = _interp(c.steps, anchors)

    tau = np.empty((c.steps, n), dtype=float)
    for i in range(n):
        phase1 = rng.uniform(0.0, 1.0)
        phase2 = rng.uniform(0.0, 1.0)
        amp1 = rng.uniform(0.08, 0.22)
        amp2 = rng.uniform(0.03, 0.10)
        local = np.exp(
            amp1 * np.sin(2.0 * np.pi * (rng.uniform(0.7, 1.7) * t + phase1))
            + amp2 * np.sin(
                2.0 * np.pi * (rng.uniform(2.0, 3.8) * t + phase2)
            )
        )
        tau[:, i] = base_tau[i] * global_scale * local

    logits = np.empty((c.steps, n), dtype=float)
    for i in range(n):
        phase = rng.uniform(0.0, 1.0)
        phase2 = rng.uniform(0.0, 1.0)
        logits[:, i] = (
            rng.uniform(0.25, 0.75)
            * np.sin(2.0 * np.pi * (rng.uniform(0.45, 1.15) * t + phase))
            + rng.uniform(0.15, 0.45)
            * np.cos(2.0 * np.pi * (rng.uniform(1.1, 2.2) * t + phase2))
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

    return {"x": x, "y": y, "horizon": horizon}


def sample_candidate_rules(config: MetaConfig) -> list[dict[str, float]]:
    """Deterministic outer-loop population of candidate developmental laws."""
    rng = np.random.default_rng(config.meta_seed)
    rows: list[dict[str, float]] = []

    # Hand rule from Gate 6 is explicitly included as a ruler.
    rows.append(
        {
            "fast_rate": 3.0,
            "slow_rate": 0.055,
            "fast_retention": 0.996,
            "readout_rate": 0.075,
        }
    )

    while len(rows) < config.candidate_count:
        rows.append(
            {
                "fast_rate": float(
                    np.exp(rng.uniform(np.log(0.35), np.log(6.0)))
                ),
                "slow_rate": float(
                    np.exp(rng.uniform(np.log(0.006), np.log(0.16)))
                ),
                "fast_retention": float(rng.uniform(0.985, 0.999)),
                "readout_rate": float(
                    np.exp(rng.uniform(np.log(0.025), np.log(0.22)))
                ),
            }
        )
    return rows


def _rich_config(meta: MetaConfig, theta: dict[str, float]) -> RichConfig:
    return RichConfig(
        steps=meta.steps,
        burnin=meta.burnin,
        seeds=1,
        moving_sizes=(meta.modes,),
        context_windows=(meta.context_window,),
        fast_rate=theta["fast_rate"],
        slow_rate=theta["slow_rate"],
        fast_retention=theta["fast_retention"],
        readout_rate=theta["readout_rate"],
        context_rate=0.85,
    )


def run_rule_on_world(
    theta: dict[str, float],
    world: dict[str, np.ndarray],
    config: MetaConfig,
    *,
    move_operator: bool = True,
    keep_trace: bool = False,
) -> dict[str, object]:
    rc = _rich_config(config, theta)
    model = BudgetMovingOperator(
        config.modes, rc, move_operator=move_operator
    )
    x = world["x"]
    y = world["y"]
    pred = np.zeros(config.steps, dtype=float)
    tau = np.zeros(config.steps, dtype=float)

    for t in range(config.steps):
        p = model.predict(float(x[t]))
        pred[t] = p
        model.learn(float(y[t]), p)
        tau[t] = model.effective_tau()

    active = slice(config.burnin, None)
    mse = float(np.mean((pred[active] - y[active]) ** 2))
    corr = float(
        np.corrcoef(tau[active], world["horizon"][active])[0, 1]
    )

    result: dict[str, object] = {
        "mse": mse,
        "horizon_correlation": corr,
        "final_effective_tau": float(tau[-1]),
    }
    if keep_trace:
        stride = 50
        result["trace"] = {
            "step": list(range(0, config.steps, stride)),
            "teacher_horizon": world["horizon"][::stride].tolist(),
            "operator_horizon": tau[::stride].tolist(),
        }
    return result


def run_context_on_world(
    world: dict[str, np.ndarray],
    config: MetaConfig,
) -> float:
    model = ContextLMS(config.context_window, rate=0.85)
    pred = np.zeros(config.steps, dtype=float)
    for t in range(config.steps):
        p, v = model.predict(float(world["x"][t]))
        pred[t] = p
        model.learn(float(world["y"][t]), p, v)
    active = slice(config.burnin, None)
    return float(
        np.mean((pred[active] - world["y"][active]) ** 2)
    )


def _mean(values) -> float:
    return float(np.mean(np.asarray(tuple(values), dtype=float)))


def _std(values) -> float:
    return float(np.std(np.asarray(tuple(values), dtype=float)))


def run_meta_battery(config: MetaConfig | None = None) -> dict:
    c = config or MetaConfig()

    train_seeds = [100 + i for i in range(c.train_worlds)]
    test_seeds = [1000 + i for i in range(c.test_worlds)]
    train_world = {s: meta_world(s, c) for s in train_seeds}
    test_world = {s: meta_world(s, c) for s in test_seeds}

    candidates = sample_candidate_rules(c)

    train_scores = []
    for theta in candidates:
        score = _mean(
            run_rule_on_world(theta, train_world[s], c)["mse"]
            for s in train_seeds
        )
        train_scores.append(score)

    selected_idx = int(np.argmin(train_scores))
    selected = candidates[selected_idx]
    hand = candidates[0]

    candidate_test_scores = []
    candidate_test_corrs = []
    for theta in candidates:
        rows = [
            run_rule_on_world(theta, test_world[s], c)
            for s in test_seeds
        ]
        candidate_test_scores.append(_mean(row["mse"] for row in rows))
        candidate_test_corrs.append(
            _mean(row["horizon_correlation"] for row in rows)
        )

    selected_rows = [
        run_rule_on_world(
            selected,
            test_world[s],
            c,
            keep_trace=(i == 0),
        )
        for i, s in enumerate(test_seeds)
    ]
    hand_rows = [
        run_rule_on_world(hand, test_world[s], c)
        for s in test_seeds
    ]
    frozen_rows = [
        run_rule_on_world(
            selected, test_world[s], c, move_operator=False
        )
        for s in test_seeds
    ]

    # Cheating upper bound: choose theta separately using each test world's
    # own labels. This is not a deployable model; it is a distance-to-oracle
    # ruler for the one-rule-for-many-worlds claim.
    oracle_scores = []
    for s in test_seeds:
        scores = [
            run_rule_on_world(theta, test_world[s], c)["mse"]
            for theta in candidates
        ]
        oracle_scores.append(float(np.min(scores)))

    context_scores = [
        run_context_on_world(test_world[s], c) for s in test_seeds
    ]

    selected_mse = _mean(row["mse"] for row in selected_rows)
    hand_mse = _mean(row["mse"] for row in hand_rows)
    frozen_mse = _mean(row["mse"] for row in frozen_rows)
    oracle_mse = _mean(oracle_scores)
    context_mse = _mean(context_scores)

    test_rank = int(
        np.argsort(np.asarray(candidate_test_scores)).tolist().index(
            selected_idx
        )
        + 1
    )

    return {
        "config": c.as_dict(),
        "train_world_seeds": train_seeds,
        "test_world_seeds": test_seeds,
        "selected_candidate_index": selected_idx,
        "selected_theta": selected,
        "hand_theta": hand,
        "selected_train_mse": float(train_scores[selected_idx]),
        "selected_test_mse": selected_mse,
        "selected_test_std": _std(row["mse"] for row in selected_rows),
        "selected_test_horizon_correlation": _mean(
            row["horizon_correlation"] for row in selected_rows
        ),
        "hand_test_mse": hand_mse,
        "frozen_test_mse": frozen_mse,
        "median_candidate_test_mse": float(
            np.median(candidate_test_scores)
        ),
        "per_world_oracle_test_mse": oracle_mse,
        "context64_test_mse": context_mse,
        "selected_to_oracle_ratio":
            selected_mse / (oracle_mse + 1e-12),
        "selected_to_hand_ratio":
            selected_mse / (hand_mse + 1e-12),
        "selected_to_median_candidate_ratio":
            selected_mse
            / (float(np.median(candidate_test_scores)) + 1e-12),
        "selected_candidate_test_rank":
            test_rank,
        "candidate_count": len(candidates),
        "test_final_operator_tau_std": float(
            np.std(
                [row["final_effective_tau"] for row in selected_rows]
            )
        ),
        "seed0_trace": selected_rows[0]["trace"],
    }


def write_json(path: str | Path, payload: dict) -> None:
    import json

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
