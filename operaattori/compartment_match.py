"""Matched-site utilities for Gate 15 branch-compartment audit."""
from __future__ import annotations

import numpy as np


def branch_runs(parents: np.ndarray, edge_active: np.ndarray, clamped: np.ndarray) -> list[list[int]]:
    """Maximal unbranched runs between clamp roots, bifurcations, and tips."""
    parents = np.asarray(parents, dtype=np.int64)
    edge_active = np.asarray(edge_active, dtype=bool)
    clamped = np.asarray(clamped, dtype=bool)

    children: dict[int, list[int]] = {}
    for i in range(1, len(parents)):
        if edge_active[i]:
            children.setdefault(int(parents[i]), []).append(i)

    starts = set(int(i) for i in np.flatnonzero(clamped))
    starts.update(
        i for i in children
        if not clamped[i] and len(children.get(i, [])) != 1
    )

    runs: list[list[int]] = []
    for start in sorted(starts):
        for child in children.get(start, []):
            path = [int(child)]
            cur = int(child)
            while len(children.get(cur, [])) == 1:
                cur = int(children[cur][0])
                path.append(cur)
            runs.append(path)
    return runs


def run_id_map(n_nodes: int, runs: list[list[int]]) -> np.ndarray:
    out = np.full(n_nodes, -1, dtype=np.int64)
    for rid, run in enumerate(runs):
        out[np.asarray(run, dtype=int)] = rid
    return out


def select_even_sites(run: list[int] | np.ndarray, count: int) -> np.ndarray:
    nodes = np.asarray(run, dtype=int)
    count = min(int(count), len(nodes))
    if count <= 0:
        return np.empty(0, dtype=int)
    idx = np.unique(np.rint(np.linspace(0, len(nodes) - 1, count)).astype(int))
    return nodes[idx]


def _features(
    sites: np.ndarray,
    zdrive_mohm: np.ndarray,
    soma_transfer: np.ndarray,
) -> np.ndarray:
    sites = np.asarray(sites, dtype=int)
    z = np.maximum(np.abs(np.asarray(zdrive_mohm)[sites]), 1e-12)
    t = np.maximum(np.abs(np.asarray(soma_transfer)[sites]), 1e-12)
    return np.column_stack([np.log(z), np.log(t)])


def greedy_dispersed_match(
    target_sites: np.ndarray,
    target_run_id: int,
    pool_sites: np.ndarray,
    run_ids: np.ndarray,
    zdrive_mohm: np.ndarray,
    soma_transfer: np.ndarray,
    *,
    min_distinct_runs: int = 4,
    reuse_penalty: float = 0.35,
) -> tuple[np.ndarray, dict]:
    """Match each target site to a passive-similar site on other branch runs.

    Matching coordinates are log driving-point impedance and log absolute soma
    transfer. Sites are unique. Early matches preferentially open distinct
    branch runs; later matches may reuse them with a mild occupancy penalty.
    """
    target_sites = np.asarray(target_sites, dtype=int)
    pool_sites = np.asarray(pool_sites, dtype=int)
    run_ids = np.asarray(run_ids, dtype=np.int64)

    valid = (
        (run_ids[pool_sites] >= 0)
        & (run_ids[pool_sites] != int(target_run_id))
        & (np.abs(np.asarray(zdrive_mohm)[pool_sites]) > 0)
        & (np.abs(np.asarray(soma_transfer)[pool_sites]) > 0)
    )
    pool = pool_sites[valid]
    if len(pool) < len(target_sites):
        raise ValueError("insufficient dispersed matching pool")

    pool_feat = _features(pool, zdrive_mohm, soma_transfer)
    target_feat = _features(target_sites, zdrive_mohm, soma_transfer)
    scale = np.std(pool_feat, axis=0)
    scale = np.where(scale > 1e-8, scale, 1.0)

    available = np.ones(len(pool), dtype=bool)
    use_by_run: dict[int, int] = {}
    chosen: list[int] = []

    for tf in target_feat:
        candidate_idx = np.flatnonzero(available)
        if not len(candidate_idx):
            raise ValueError("matching pool exhausted")

        # Until enough branch runs have been opened, prefer an unused run if
        # one is available. This prevents "dispersed" from becoming another
        # single-branch cluster merely because that branch matches best.
        used_runs = set(use_by_run)
        if len(used_runs) < min_distinct_runs:
            unused = np.asarray(
                [run_ids[pool[k]] not in used_runs for k in candidate_idx],
                dtype=bool,
            )
            if np.any(unused):
                candidate_idx = candidate_idx[unused]

        delta = (pool_feat[candidate_idx] - tf) / scale
        score = np.sum(delta * delta, axis=1)
        score += np.asarray(
            [
                reuse_penalty * use_by_run.get(int(run_ids[pool[k]]), 0)
                for k in candidate_idx
            ],
            dtype=float,
        )
        best_local = int(np.argmin(score))
        k = int(candidate_idx[best_local])
        site = int(pool[k])
        chosen.append(site)
        available[k] = False
        rid = int(run_ids[site])
        use_by_run[rid] = use_by_run.get(rid, 0) + 1

    chosen_arr = np.asarray(chosen, dtype=int)
    target_f = _features(target_sites, zdrive_mohm, soma_transfer)
    match_f = _features(chosen_arr, zdrive_mohm, soma_transfer)
    log_abs_error = np.abs(match_f - target_f)

    diagnostics = {
        "distinct_match_runs": int(len(set(int(run_ids[s]) for s in chosen_arr))),
        "median_abs_log_z_error": float(np.median(log_abs_error[:, 0])),
        "median_abs_log_transfer_error": float(np.median(log_abs_error[:, 1])),
        "median_z_ratio_factor": float(np.exp(np.median(log_abs_error[:, 0]))),
        "median_transfer_ratio_factor": float(np.exp(np.median(log_abs_error[:, 1]))),
    }
    return chosen_arr, diagnostics


def normalized_offdiagonal_coupling(z_mohm: np.ndarray) -> float:
    Z = np.asarray(z_mohm, dtype=float)
    if len(Z) < 2:
        return 0.0
    d = np.maximum(np.diag(Z), 1e-12)
    denom = np.sqrt(d[:, None] * d[None, :])
    C = np.abs(Z) / denom
    mask = ~np.eye(len(Z), dtype=bool)
    return float(np.median(C[mask]))
