from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

FCI_COMMIT = "75ad8b4d81a7f51bf888b30650c543592340db06"

DT_MS = 0.05
IMPULSE_MS = 5.0
POST_MS = 60.0
TSTOP_MS = IMPULSE_MS + POST_MS
IMPULSE_NA = 0.001
N_BRANCHES = 6
SITE_X = (0.25, 0.50, 0.75)
PCA_COMPONENTS = 8
RIDGE_ALPHA = 1.0
V_REST_MV = -70.0

COMMON_MODEL = (
    "simulating_neurons/neuron_models/rat/hay/"
    "Rat_L5b_PC_2_Hay_passive_dends_simple_soma"
)
HOC_MODEL = "simulating_neurons/neuron_models/passive_dends_simple_soma_model.hoc"


def git_head(path: Path) -> str:
    import subprocess

    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        text=True,
    ).strip()


def load_panel(fci_root: Path) -> list[dict]:
    rows = []
    with (fci_root / "models.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "order": int(row["order_in_fig2"]),
                    "species": row["species"],
                    "layer": row["layer"],
                    "morphology_identifier": row["morphology_identifier"],
                    "model_folder": row["model_folder"],
                    "morphology_file": row["morphology_file"],
                }
            )
    rows.sort(key=lambda x: x["order"])
    if len(rows) != 24:
        raise RuntimeError(f"expected 24 FCI models, got {len(rows)}")
    return rows


def setup_neuron(fci_root: Path) -> None:
    import neuron
    from neuron import h

    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")
    neuron.load_mechanisms(str(fci_root / COMMON_MODEL))
    h.load_file(str(fci_root / HOC_MODEL))
    h.celsius = 34.0


def instantiate_matched_passive(
    fci_root: Path,
    model: dict,
):
    from neuron import h

    model_dir = fci_root / model["model_folder"]
    morph_dir = model_dir / "morphologies"
    morph_file = model["morphology_file"]

    cell = h.PassiveDendsSimpleSomaModel(
        str(morph_dir),
        str(morph_file),
        25.0,
    )

    # Remove fitted active soma/axon differences and impose one electrical
    # regime so this is morphology -> passive operator only.
    for sec in cell.all:
        for mech in ("kv", "na"):
            try:
                sec.uninsert(mech)
            except Exception:
                pass
        sec.insert("pas")
        sec.Ra = 150.0
        sec.cm = 1.0
        sec.g_pas = 1.0 / 20000.0
        sec.e_pas = V_REST_MV

    h.distance(0, 0.5, sec=cell.soma[0])
    return cell


def section_name(sec) -> str:
    return str(sec.name())


def subtree_sections(root_sec) -> list:
    out = []
    stack = [root_sec]
    seen = set()
    while stack:
        sec = stack.pop()
        name = section_name(sec)
        if name in seen:
            continue
        seen.add(name)
        out.append(sec)
        try:
            children = list(sec.children())
        except Exception:
            children = []
        stack.extend(children)
    return out


def branch_order(sec) -> int:
    order = 0
    current = sec
    seen = set()
    while True:
        name = section_name(current)
        if name in seen:
            break
        seen.add(name)
        try:
            parent_seg = current.parentseg()
        except Exception:
            parent_seg = None
        if parent_seg is None:
            break
        parent_sec = parent_seg.sec
        if "soma" in section_name(parent_sec):
            order += 1
            break
        order += 1
        current = parent_sec
        if order > 1000:
            raise RuntimeError("branch-order traversal did not terminate")
    return order


def section_area_um2(sec) -> float:
    from neuron import h

    total = 0.0
    for seg in sec:
        total += float(h.area(float(seg.x), sec=sec))
    return total


def safe_log(value: float, floor: float = 1e-9) -> float:
    return float(math.log(max(float(value), floor)))


def cell_global_geometry(cell) -> dict:
    from neuron import h

    apical = list(cell.apical)
    basal = list(cell.basal)
    dendritic = apical + basal
    if not dendritic:
        raise RuntimeError("cell has no dendritic sections")

    total_length = float(sum(float(sec.L) for sec in dendritic))
    apical_length = float(sum(float(sec.L) for sec in apical))
    total_area = float(sum(section_area_um2(sec) for sec in dendritic))

    max_path = 0.0
    for sec in dendritic:
        for x in (0.0, 0.5, 1.0):
            try:
                d = float(h.distance(float(x), sec=sec))
            except Exception:
                d = 0.0
            if np.isfinite(d):
                max_path = max(max_path, d)

    return {
        "total_dendritic_length_um": total_length,
        "total_dendritic_area_um2": total_area,
        "max_dendritic_path_um": max_path,
        "total_apical_length_um": apical_length,
    }


def choose_branches(cell) -> list:
    apical = list(cell.apical)
    if len(apical) < N_BRANCHES:
        raise RuntimeError(
            f"need {N_BRANCHES} apical sections, found {len(apical)}"
        )
    apical.sort(
        key=lambda sec: (-float(sec.L), section_name(sec))
    )
    return apical[:N_BRANCHES]


def branch_record(cell, sec, branch_index: int, global_geom: dict) -> dict:
    from neuron import h

    subtree = subtree_sections(sec)
    subtree_length = float(sum(float(s.L) for s in subtree))
    diameters = [float(sec(x).diam) for x in SITE_X]
    midpoint_path = float(h.distance(0.5, sec=sec))

    local_features = [
        safe_log(float(sec.L)),
        safe_log(float(np.mean(diameters))),
        safe_log(midpoint_path),
        float(branch_order(sec)),
        safe_log(subtree_length),
        safe_log(1.0 + len(subtree)),
    ]
    global_features = [
        safe_log(global_geom["total_dendritic_length_um"]),
        safe_log(global_geom["total_dendritic_area_um2"]),
        safe_log(global_geom["max_dendritic_path_um"]),
        safe_log(global_geom["total_apical_length_um"]),
    ]

    return {
        "branch_index": int(branch_index),
        "section": section_name(sec),
        "sec": sec,
        "site_x": list(SITE_X),
        "local_features": local_features,
        "global_features": global_features,
        "features": local_features + global_features,
        "descriptors": {
            "section_length_um": float(sec.L),
            "mean_diameter_um": float(np.mean(diameters)),
            "midpoint_path_um": midpoint_path,
            "branch_order": int(branch_order(sec)),
            "subtree_length_um": subtree_length,
            "subtree_sections": int(len(subtree)),
            **global_geom,
        },
    }


def run_noinput(cell, branch_rows: list[dict]) -> dict:
    from neuron import h

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
    local_vecs = []
    local_keys = []
    for branch in branch_rows:
        sec = branch["sec"]
        for site_index, x in enumerate(SITE_X):
            local_vecs.append(h.Vector().record(sec(x)._ref_v))
            local_keys.append((branch["branch_index"], site_index))

    h.dt = DT_MS
    h.finitialize(V_REST_MV)
    h.fcurrent()
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    post = (
        (t >= IMPULSE_MS)
        & (t <= IMPULSE_MS + POST_MS + 1e-9)
    )
    tp = t[post] - IMPULSE_MS

    local = {}
    for key, vec in zip(local_keys, local_vecs):
        local[key] = np.asarray(vec, dtype=float)[post]

    return {
        "t": tp,
        "soma": np.asarray(soma_vec, dtype=float)[post],
        "local": local,
    }


def measure_branch_operator(
    cell,
    branch: dict,
    noinput: dict,
) -> tuple[np.ndarray, np.ndarray]:
    from neuron import h

    sec = branch["sec"]
    n = len(noinput["t"])
    G = np.zeros((3, 3, n), dtype=float)
    T = np.zeros((3, n), dtype=float)

    for source_index, source_x in enumerate(SITE_X):
        stim = h.IClamp(float(source_x), sec=sec)
        stim.delay = IMPULSE_MS
        stim.dur = DT_MS
        stim.amp = IMPULSE_NA

        tvec = h.Vector().record(h._ref_t)
        soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
        local_vecs = [
            h.Vector().record(sec(x)._ref_v)
            for x in SITE_X
        ]

        h.dt = DT_MS
        h.finitialize(V_REST_MV)
        h.fcurrent()
        h.continuerun(TSTOP_MS)

        t = np.asarray(tvec, dtype=float)
        post = (
            (t >= IMPULSE_MS)
            & (t <= IMPULSE_MS + POST_MS + 1e-9)
        )
        tp = t[post] - IMPULSE_MS
        if not np.allclose(
            tp, noinput["t"], rtol=0, atol=1e-12
        ):
            raise RuntimeError("operator time grid differs from no-input grid")

        soma = np.asarray(soma_vec, dtype=float)[post]
        T[source_index] = (
            soma - noinput["soma"]
        ) / IMPULSE_NA

        for target_index, vec in enumerate(local_vecs):
            local = np.asarray(vec, dtype=float)[post]
            base = noinput["local"][
                (branch["branch_index"], target_index)
            ]
            G[target_index, source_index] = (
                local - base
            ) / IMPULSE_NA

        stim.amp = 0.0

    if not np.all(np.isfinite(G)) or not np.all(np.isfinite(T)):
        raise FloatingPointError("non-finite cross-cell operator")
    return G, T


def collect_panel(fci_root: Path, panel: list[dict]) -> dict:
    rows = []
    time_grid = None
    cell_receipts = []

    for ci, model in enumerate(panel):
        cell = instantiate_matched_passive(fci_root, model)
        global_geom = cell_global_geometry(cell)
        secs = choose_branches(cell)
        branches = [
            branch_record(cell, sec, bi, global_geom)
            for bi, sec in enumerate(secs)
        ]
        noinput = run_noinput(cell, branches)

        if time_grid is None:
            time_grid = noinput["t"].copy()
        elif not np.allclose(
            time_grid, noinput["t"], rtol=0, atol=1e-12
        ):
            raise RuntimeError("cell time grids differ")

        cell_receipt = {
            "cell_order": int(model["order"]),
            "species": model["species"],
            "layer": model["layer"],
            "morphology_identifier": model["morphology_identifier"],
            "model_folder": model["model_folder"],
            "morphology_file": model["morphology_file"],
            "branches": [],
        }

        for branch in branches:
            G, T = measure_branch_operator(
                cell, branch, noinput
            )
            row = {
                "cell_order": int(model["order"]),
                "species": model["species"],
                "layer": model["layer"],
                "morphology_identifier": model["morphology_identifier"],
                "branch_index": int(branch["branch_index"]),
                "section": branch["section"],
                "features": np.asarray(
                    branch["features"], dtype=float
                ),
                "local_features": np.asarray(
                    branch["local_features"], dtype=float
                ),
                "descriptors": branch["descriptors"],
                "G": G.reshape(-1),
                "T": T.reshape(-1),
            }
            rows.append(row)
            cell_receipt["branches"].append(
                {
                    "branch_index": int(branch["branch_index"]),
                    "section": branch["section"],
                    "descriptors": branch["descriptors"],
                }
            )

        cell_receipts.append(cell_receipt)
        print(
            f"[{ci+1:02d}/24] {model['species']:5s} "
            f"{model['morphology_identifier']:>12s} "
            f"branches={len(branches)} "
            f"apical_maxL={max(float(s.L) for s in secs):.1f}um"
        )

    if len(rows) != 24 * N_BRANCHES:
        raise RuntimeError(
            f"expected {24*N_BRANCHES} branch operators, got {len(rows)}"
        )

    return {
        "rows": rows,
        "time": time_grid,
        "cells": cell_receipts,
    }


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def nrmse(actual: np.ndarray, pred: np.ndarray) -> float:
    return rms(np.asarray(pred) - np.asarray(actual)) / (
        rms(actual) + 1e-30
    )


def standardize(
    X_train: np.ndarray,
    X_test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(X_train, axis=0)
    std = np.std(X_train, axis=0)
    std = np.where(std < 1e-12, 1.0, std)
    return (
        (X_train - mean) / std,
        (X_test - mean) / std,
    )


def pca_basis(Y_train: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(Y_train, axis=0)
    Yc = Y_train - mean

    # Economy PCA through the sample Gram matrix.  This avoids an expensive
    # full SVD of the long time-domain operator vector in every LOCO fold.
    gram = Yc @ Yc.T
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]

    basis_rows = []
    for idx in range(min(k, len(values))):
        value = float(values[idx])
        if value <= 1e-18:
            continue
        vec = vectors[:, idx]
        basis = (vec @ Yc) / math.sqrt(value)
        norm = np.linalg.norm(basis)
        if norm <= 1e-15:
            continue
        basis_rows.append(basis / norm)

    if not basis_rows:
        raise RuntimeError("PCA basis is empty")
    return mean, np.stack(basis_rows, axis=0)


def ridge_scores(
    X_train_z: np.ndarray,
    scores_train: np.ndarray,
    X_test_z: np.ndarray,
) -> np.ndarray:
    X1 = np.column_stack(
        [np.ones(len(X_train_z)), X_train_z]
    )
    Xt = np.column_stack(
        [np.ones(len(X_test_z)), X_test_z]
    )
    reg = np.eye(X1.shape[1], dtype=float) * RIDGE_ALPHA
    reg[0, 0] = 0.0
    beta = np.linalg.solve(
        X1.T @ X1 + reg,
        X1.T @ scores_train,
    )
    return Xt @ beta


def fit_predict_operator(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    X_test: np.ndarray,
    Y_test_actual: np.ndarray,
) -> dict:
    X_train_z, X_test_z = standardize(X_train, X_test)
    mean, basis = pca_basis(
        Y_train, PCA_COMPONENTS
    )
    train_scores = (Y_train - mean) @ basis.T
    pred_scores = ridge_scores(
        X_train_z,
        train_scores,
        X_test_z,
    )
    predicted = mean + pred_scores @ basis

    oracle_scores = (Y_test_actual - mean) @ basis.T
    pca_oracle = mean + oracle_scores @ basis

    # Baselines are generated from training data only.
    mean_pred = np.repeat(
        mean[None, :], len(X_test), axis=0
    )

    nearest_idx = []
    for row in X_test_z:
        distances = np.sum(
            (X_train_z - row[None, :]) ** 2,
            axis=1,
        )
        nearest_idx.append(int(np.argmin(distances)))
    nearest_pred = Y_train[
        np.asarray(nearest_idx, dtype=int)
    ]

    return {
        "predicted": predicted,
        "pca_oracle": pca_oracle,
        "mean_pred": mean_pred,
        "nearest_pred": nearest_pred,
        "basis_rank": int(basis.shape[0]),
    }


def row_matrix(rows: list[dict], key: str) -> np.ndarray:
    return np.stack(
        [np.asarray(row[key], dtype=float) for row in rows],
        axis=0,
    )


def score_predictions(
    rows_test: list[dict],
    predG: dict,
    predT: dict,
    model_name: str,
) -> list[dict]:
    actual_G = row_matrix(rows_test, "G")
    actual_T = row_matrix(rows_test, "T")

    if model_name == "morphology":
        Gp = predG["predicted"]
        Tp = predT["predicted"]
    elif model_name == "mean":
        Gp = predG["mean_pred"]
        Tp = predT["mean_pred"]
    elif model_name == "nearest":
        Gp = predG["nearest_pred"]
        Tp = predT["nearest_pred"]
    elif model_name == "pca_oracle":
        Gp = predG["pca_oracle"]
        Tp = predT["pca_oracle"]
    else:
        raise ValueError(model_name)

    scored = []
    for i, row in enumerate(rows_test):
        eg = nrmse(actual_G[i], Gp[i])
        et = nrmse(actual_T[i], Tp[i])
        scored.append(
            {
                "cell_order": int(row["cell_order"]),
                "species": row["species"],
                "morphology_identifier": row["morphology_identifier"],
                "branch_index": int(row["branch_index"]),
                "section": row["section"],
                "G_nrmse": float(eg),
                "T_nrmse": float(et),
                "joint_nrmse": float(0.5 * (eg + et)),
            }
        )
    return scored


def aggregate_scores(scores: list[dict]) -> dict:
    g = np.asarray([x["G_nrmse"] for x in scores], dtype=float)
    t = np.asarray([x["T_nrmse"] for x in scores], dtype=float)
    j = np.asarray([x["joint_nrmse"] for x in scores], dtype=float)

    per_cell = {}
    for cell_order in sorted({x["cell_order"] for x in scores}):
        vals = [
            x for x in scores
            if x["cell_order"] == cell_order
        ]
        per_cell[str(cell_order)] = {
            "species": vals[0]["species"],
            "morphology_identifier": vals[0]["morphology_identifier"],
            "median_G_nrmse": float(
                np.median([x["G_nrmse"] for x in vals])
            ),
            "median_T_nrmse": float(
                np.median([x["T_nrmse"] for x in vals])
            ),
            "median_joint_nrmse": float(
                np.median([x["joint_nrmse"] for x in vals])
            ),
        }

    return {
        "branches": int(len(scores)),
        "median_G_nrmse": float(np.median(g)),
        "median_T_nrmse": float(np.median(t)),
        "median_joint_nrmse": float(np.median(j)),
        "per_cell": per_cell,
    }


def loco(panel_rows: list[dict], feature_key: str) -> dict:
    all_scores = {
        "morphology": [],
        "mean": [],
        "nearest": [],
        "pca_oracle": [],
    }
    fold_receipts = []

    cell_orders = sorted(
        {int(row["cell_order"]) for row in panel_rows}
    )

    for fold_i, heldout in enumerate(cell_orders):
        train = [
            row for row in panel_rows
            if int(row["cell_order"]) != heldout
        ]
        test = [
            row for row in panel_rows
            if int(row["cell_order"]) == heldout
        ]

        X_train = row_matrix(train, feature_key)
        X_test = row_matrix(test, feature_key)
        YG_train = row_matrix(train, "G")
        YG_test = row_matrix(test, "G")
        YT_train = row_matrix(train, "T")
        YT_test = row_matrix(test, "T")

        predG = fit_predict_operator(
            X_train, YG_train, X_test, YG_test
        )
        predT = fit_predict_operator(
            X_train, YT_train, X_test, YT_test
        )

        fold_result = {
            "heldout_cell_order": int(heldout),
            "species": test[0]["species"],
            "morphology_identifier": test[0][
                "morphology_identifier"
            ],
            "G_pca_rank": predG["basis_rank"],
            "T_pca_rank": predT["basis_rank"],
        }

        for name in all_scores:
            scores = score_predictions(
                test, predG, predT, name
            )
            all_scores[name].extend(scores)
            fold_result[name] = {
                "median_joint_nrmse": float(
                    np.median(
                        [x["joint_nrmse"] for x in scores]
                    )
                ),
                "median_G_nrmse": float(
                    np.median([x["G_nrmse"] for x in scores])
                ),
                "median_T_nrmse": float(
                    np.median([x["T_nrmse"] for x in scores])
                ),
            }

        print(
            f"LOCO [{fold_i+1:02d}/24] "
            f"{heldout:02d} {test[0]['species']:5s} "
            f"{test[0]['morphology_identifier']:>12s} "
            f"morph={fold_result['morphology']['median_joint_nrmse']:.3f} "
            f"near={fold_result['nearest']['median_joint_nrmse']:.3f} "
            f"mean={fold_result['mean']['median_joint_nrmse']:.3f} "
            f"pca={fold_result['pca_oracle']['median_joint_nrmse']:.3f}"
        )
        fold_receipts.append(fold_result)

    aggregates = {
        name: aggregate_scores(scores)
        for name, scores in all_scores.items()
    }

    wins_nearest = 0
    for fold in fold_receipts:
        if (
            fold["morphology"]["median_joint_nrmse"]
            < fold["nearest"]["median_joint_nrmse"]
        ):
            wins_nearest += 1

    return {
        "feature_key": feature_key,
        "folds": fold_receipts,
        "aggregate": aggregates,
        "morphology_beats_nearest_cells": int(wins_nearest),
    }


def species_transfer(
    panel_rows: list[dict],
    train_species: str,
    test_species: str,
) -> dict:
    train = [
        row for row in panel_rows
        if row["species"] == train_species
    ]
    test = [
        row for row in panel_rows
        if row["species"] == test_species
    ]

    X_train = row_matrix(train, "features")
    X_test = row_matrix(test, "features")
    predG = fit_predict_operator(
        X_train,
        row_matrix(train, "G"),
        X_test,
        row_matrix(test, "G"),
    )
    predT = fit_predict_operator(
        X_train,
        row_matrix(train, "T"),
        X_test,
        row_matrix(test, "T"),
    )

    scored = {}
    for name in (
        "morphology",
        "mean",
        "nearest",
        "pca_oracle",
    ):
        scored[name] = aggregate_scores(
            score_predictions(
                test, predG, predT, name
            )
        )

    return {
        "train_species": train_species,
        "test_species": test_species,
        "train_branches": len(train),
        "test_branches": len(test),
        "scores": scored,
    }


def summarize_by_species(
    loco_result: dict,
) -> dict:
    folds = loco_result["folds"]
    out = {}
    for species in ("rat", "human"):
        vals = [
            fold["morphology"]["median_joint_nrmse"]
            for fold in folds
            if fold["species"] == species
        ]
        out[species] = {
            "cells": int(len(vals)),
            "median_cell_joint_nrmse": float(
                np.median(vals)
            ),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/cross_cell_operator.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source is not pinned to the preregistered commit")

    setup_neuron(fci_root)
    panel = load_panel(fci_root)
    measured = collect_panel(fci_root, panel)
    rows = measured["rows"]

    full = loco(rows, "features")
    local = loco(rows, "local_features")

    agg = full["aggregate"]
    morph = agg["morphology"]
    mean = agg["mean"]
    nearest = agg["nearest"]
    pca = agg["pca_oracle"]

    mean_ratio = (
        morph["median_joint_nrmse"]
        / (mean["median_joint_nrmse"] + 1e-30)
    )
    nearest_ratio = (
        morph["median_joint_nrmse"]
        / (nearest["median_joint_nrmse"] + 1e-30)
    )

    primary_pass = (
        morph["median_joint_nrmse"] <= 0.20
        and morph["median_T_nrmse"] <= 0.20
        and morph["median_G_nrmse"] <= 0.20
        and mean_ratio <= 0.80
        and nearest_ratio <= 0.90
        and full["morphology_beats_nearest_cells"] >= 16
    )

    if primary_pass:
        classification = (
            "CROSS_CELL_OPERATOR_PREDICTABLE_FROM_MORPHOLOGY"
        )
        interpretation = (
            "A fixed morphology-only linear map predicts the local Green and "
            "site-to-soma operator packs of held-out released neuron "
            "morphologies better than both mean and nearest-morphology "
            "attackers under leave-one-cell-out validation."
        )
    elif pca["median_joint_nrmse"] <= 0.10:
        classification = (
            "CROSS_CELL_OPERATOR_LOW_DIMENSIONAL_BUT_MORPHOLOGY_MAP_WEAK"
        )
        interpretation = (
            "Held-out operators lie within the locked low-dimensional training "
            "basis, but the preregistered morphology features do not predict "
            "their coordinates accurately enough."
        )
    else:
        classification = (
            "CROSS_CELL_OPERATOR_FAMILY_NOT_CAPTURED_BY_LOCKED_BASIS"
        )
        interpretation = (
            "Even the training-only eight-component PCA basis does not "
            "reconstruct held-out cells at the preregistered accuracy, so the "
            "cross-cell operator family is not captured by this fixed basis."
        )

    secondary = {
        "local_only_median_joint_nrmse": float(
            local["aggregate"]["morphology"][
                "median_joint_nrmse"
            ]
        ),
        "full_feature_median_joint_nrmse": float(
            morph["median_joint_nrmse"]
        ),
        "full_to_local_only_error_ratio": float(
            morph["median_joint_nrmse"]
            / (
                local["aggregate"]["morphology"][
                    "median_joint_nrmse"
                ]
                + 1e-30
            )
        ),
        "loco_species": summarize_by_species(full),
        "rat_to_human": species_transfer(
            rows, "rat", "human"
        ),
        "human_to_rat": species_transfer(
            rows, "human", "rat"
        ),
    }

    summary = {
        "object": (
            "leave-one-cell-out prediction of matched-passive local Green and "
            "site-to-soma transport operators from morphology-only features"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "cells": 24,
            "branches_per_cell": N_BRANCHES,
            "sites_per_branch": len(SITE_X),
            "site_x": list(SITE_X),
            "matched_passive": {
                "Ra_ohm_cm": 150.0,
                "Cm_uF_cm2": 1.0,
                "Rm_ohm_cm2": 20000.0,
                "E_pas_mV": V_REST_MV,
                "active_Na_K_removed": True,
            },
            "dt_ms": DT_MS,
            "impulse_nA": IMPULSE_NA,
            "impulse_ms": IMPULSE_MS,
            "post_ms": POST_MS,
            "pca_components": PCA_COMPONENTS,
            "ridge_alpha": RIDGE_ALPHA,
            "features": [
                "log_section_length",
                "log_mean_diameter",
                "log_midpoint_path_distance",
                "branch_order",
                "log_subtree_total_length",
                "log_1plus_subtree_section_count",
                "log_total_dendritic_length",
                "log_total_dendritic_area",
                "log_max_dendritic_path",
                "log_total_apical_length",
            ],
            "no_species_feature": True,
            "no_fci_feature": True,
            "heldout_unit": "whole cell",
            "thresholds_locked_before_run": {
                "median_joint_nrmse_max": 0.20,
                "median_T_nrmse_max": 0.20,
                "median_G_nrmse_max": 0.20,
                "morphology_to_mean_error_ratio_max": 0.80,
                "morphology_to_nearest_error_ratio_max": 0.90,
                "morphology_beats_nearest_cells_min": 16,
                "pca_oracle_joint_nrmse_lowdim_max": 0.10,
            },
        },
        "time_samples": int(len(measured["time"])),
        "cells": measured["cells"],
        "leave_one_cell_out_full_features": full,
        "leave_one_cell_out_local_only": local,
        "secondary": secondary,
        "primary_ratios": {
            "morphology_to_mean_joint_error_ratio": float(
                mean_ratio
            ),
            "morphology_to_nearest_joint_error_ratio": float(
                nearest_ratio
            ),
        },
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "If this fails, do not rescue it with neural networks, larger "
            "feature dictionaries, polynomial kernels or the nonlinear NMDA "
            "circuit before identifying the failure mode."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori 24-cell morphology -> operator audit")
    print()
    print(
        "morphology median joint NRMSE:   "
        f"{morph['median_joint_nrmse']:.4f}"
    )
    print(
        "  G median NRMSE:                "
        f"{morph['median_G_nrmse']:.4f}"
    )
    print(
        "  T median NRMSE:                "
        f"{morph['median_T_nrmse']:.4f}"
    )
    print(
        "mean attacker joint NRMSE:       "
        f"{mean['median_joint_nrmse']:.4f}"
    )
    print(
        "nearest attacker joint NRMSE:    "
        f"{nearest['median_joint_nrmse']:.4f}"
    )
    print(
        "PCA oracle joint NRMSE:           "
        f"{pca['median_joint_nrmse']:.4f}"
    )
    print(
        "morphology / mean:               "
        f"{mean_ratio:.4f}"
    )
    print(
        "morphology / nearest:            "
        f"{nearest_ratio:.4f}"
    )
    print(
        "held-out cells beating nearest:  "
        f"{full['morphology_beats_nearest_cells']} / 24"
    )
    print(
        "local-only joint NRMSE:           "
        f"{local['aggregate']['morphology']['median_joint_nrmse']:.4f}"
    )
    print(
        "full/local-only:                 "
        f"{secondary['full_to_local_only_error_ratio']:.4f}"
    )
    print(
        "rat LOCO median cell joint:       "
        f"{secondary['loco_species']['rat']['median_cell_joint_nrmse']:.4f}"
    )
    print(
        "human LOCO median cell joint:     "
        f"{secondary['loco_species']['human']['median_cell_joint_nrmse']:.4f}"
    )
    print(
        "rat -> human joint:              "
        f"{secondary['rat_to_human']['scores']['morphology']['median_joint_nrmse']:.4f}"
    )
    print(
        "human -> rat joint:              "
        f"{secondary['human_to_rat']['scores']['morphology']['median_joint_nrmse']:.4f}"
    )
    print(f"classification: {classification}")

    assert len(rows) == 144
    assert np.all(
        np.isfinite(
            [x["joint_nrmse"] for x in full["aggregate"]["morphology"]["per_cell"].values()]
        )
    )


if __name__ == "__main__":
    main()
