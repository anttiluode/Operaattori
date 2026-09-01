from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "audits"
if str(AUDITS) not in sys.path:
    sys.path.insert(0, str(AUDITS))

from cross_cell_operator import (
    FCI_COMMIT,
    N_BRANCHES,
    PCA_COMPONENTS,
    RIDGE_ALPHA,
    SITE_X,
    aggregate_scores,
    branch_order,
    cell_global_geometry,
    choose_branches,
    git_head,
    instantiate_matched_passive,
    load_panel,
    measure_branch_operator,
    pca_basis,
    ridge_scores,
    row_matrix,
    run_noinput,
    score_predictions,
    section_name,
    setup_neuron,
    standardize,
    subtree_sections,
)

RA_OHM_CM = 150.0
RM_OHM_CM2 = 20000.0

ELECTRO_NAMES = [
    "section_electrotonic_length",
    "soma_to_mid_electrotonic_path",
    "log_mid_characteristic_admittance_proxy",
    "log_proximal_to_distal_taper",
    "distal_rall_load_ratio",
    "proximal_rall_mismatch",
    "subtree_total_electrotonic_length",
    "subtree_max_electrotonic_depth",
    "log_terminal_admittance_mass",
    "log_rall_equivalent_terminal_diameter",
    "cell_total_dendritic_electrotonic_length",
    "cell_max_dendritic_electrotonic_path",
]


def safe_log(value: float, floor: float = 1e-12) -> float:
    return float(math.log(max(float(value), floor)))


def cable_lambda_um(diam_um: float) -> float:
    d_cm = max(float(diam_um), 1e-9) * 1e-4
    return float(
        1e4 * math.sqrt(
            d_cm * RM_OHM_CM2 / (4.0 * RA_OHM_CM)
        )
    )


def section_electrotonic_length(sec) -> float:
    total = 0.0
    nseg = max(int(sec.nseg), 1)
    ds = float(sec.L) / nseg
    for seg in sec:
        total += ds / cable_lambda_um(float(seg.diam))
    return float(total)


def proximal_diam(sec) -> float:
    return float(sec(1e-6).diam)


def distal_diam(sec) -> float:
    return float(sec(1.0 - 1e-6).diam)


def midpoint_diam(sec) -> float:
    return float(sec(0.5).diam)


def parent_section(sec):
    try:
        pseg = sec.parentseg()
    except Exception:
        pseg = None
    return None if pseg is None else pseg.sec


def parent_electrotonic_path_to_prox(sec) -> float:
    total = 0.0
    current = sec
    seen = set()
    while True:
        parent = parent_section(current)
        if parent is None:
            break
        name = section_name(parent)
        if name in seen:
            raise RuntimeError("cycle in parent traversal")
        seen.add(name)
        if "soma" in name:
            break
        total += section_electrotonic_length(parent)
        current = parent
        if len(seen) > 1000:
            raise RuntimeError("parent traversal too deep")
    return float(total)


def subtree_electrotonic_depth(root) -> float:
    memo = {}

    def rec(sec):
        key = section_name(sec)
        if key in memo:
            return memo[key]
        own = section_electrotonic_length(sec)
        children = list(sec.children())
        if not children:
            value = own
        else:
            value = own + max(rec(child) for child in children)
        memo[key] = float(value)
        return float(value)

    return rec(root)


def terminal_sections(root) -> list:
    terminals = []
    for sec in subtree_sections(root):
        if len(list(sec.children())) == 0:
            terminals.append(sec)
    return terminals


def distal_rall_load_ratio(sec) -> float:
    children = list(sec.children())
    if not children:
        return 0.0
    denom = max(distal_diam(sec), 1e-9) ** 1.5
    numer = sum(
        max(proximal_diam(child), 1e-9) ** 1.5
        for child in children
    )
    return float(numer / denom)


def proximal_rall_mismatch(sec) -> float:
    parent = parent_section(sec)
    if parent is None or "soma" in section_name(parent):
        return 1.0
    siblings = [
        child for child in list(parent.children())
        if section_name(child) != section_name(sec)
    ]
    parent_y = max(distal_diam(parent), 1e-9) ** 1.5
    child_y = max(proximal_diam(sec), 1e-9) ** 1.5
    sibling_y = sum(
        max(proximal_diam(sib), 1e-9) ** 1.5
        for sib in siblings
    )
    return float(parent_y / max(child_y + sibling_y, 1e-12))


def cell_electrotonic_context(cell) -> dict:
    dendritic = list(cell.apical) + list(cell.basal)
    sec_e = {
        section_name(sec): section_electrotonic_length(sec)
        for sec in dendritic
    }
    total = float(sum(sec_e.values()))

    max_path = 0.0
    for sec in dendritic:
        path = (
            parent_electrotonic_path_to_prox(sec)
            + sec_e[section_name(sec)]
        )
        max_path = max(max_path, path)

    return {
        "total_dendritic_electrotonic_length": total,
        "max_dendritic_electrotonic_path": float(max_path),
    }


def electrotonic_features(cell, sec, cell_ctx: dict) -> list[float]:
    subtree = subtree_sections(sec)
    subtree_total_e = float(
        sum(section_electrotonic_length(s) for s in subtree)
    )
    terminals = terminal_sections(sec)
    terminal_mass = float(
        sum(
            max(distal_diam(term), 1e-9) ** 1.5
            for term in terminals
        )
    )
    rall_eq_diam = terminal_mass ** (2.0 / 3.0)

    section_e = section_electrotonic_length(sec)
    soma_mid_e = (
        parent_electrotonic_path_to_prox(sec)
        + 0.5 * section_e
    )
    dmid = max(midpoint_diam(sec), 1e-9)
    taper = max(proximal_diam(sec), 1e-9) / max(
        distal_diam(sec), 1e-9
    )

    return [
        float(section_e),
        float(soma_mid_e),
        safe_log(dmid ** 1.5),
        float(math.log(taper)),
        float(distal_rall_load_ratio(sec)),
        float(proximal_rall_mismatch(sec)),
        float(subtree_total_e),
        float(subtree_electrotonic_depth(sec)),
        safe_log(terminal_mass),
        safe_log(rall_eq_diam),
        float(cell_ctx["total_dendritic_electrotonic_length"]),
        float(cell_ctx["max_dendritic_electrotonic_path"]),
    ]


def gross_features(sec, global_geom: dict) -> list[float]:
    subtree = subtree_sections(sec)
    subtree_length = float(sum(float(s.L) for s in subtree))
    diameters = [float(sec(x).diam) for x in SITE_X]

    from neuron import h
    midpoint_path = float(h.distance(0.5, sec=sec))

    local = [
        safe_log(float(sec.L)),
        safe_log(float(np.mean(diameters))),
        safe_log(midpoint_path),
        float(branch_order(sec)),
        safe_log(subtree_length),
        safe_log(1.0 + len(subtree)),
    ]
    global_f = [
        safe_log(global_geom["total_dendritic_length_um"]),
        safe_log(global_geom["total_dendritic_area_um2"]),
        safe_log(global_geom["max_dendritic_path_um"]),
        safe_log(global_geom["total_apical_length_um"]),
    ]
    return local + global_f


def collect_panel(fci_root: Path, panel: list[dict]) -> dict:
    rows = []
    time_grid = None
    cells = []

    for ci, model in enumerate(panel):
        cell = instantiate_matched_passive(fci_root, model)
        global_geom = cell_global_geometry(cell)
        cell_ctx = cell_electrotonic_context(cell)
        secs = choose_branches(cell)

        branch_rows = []
        for bi, sec in enumerate(secs):
            ef = electrotonic_features(cell, sec, cell_ctx)
            gf = gross_features(sec, global_geom)
            branch_rows.append(
                {
                    "branch_index": int(bi),
                    "section": section_name(sec),
                    "sec": sec,
                    "electro_features": ef,
                    "gross_features": gf,
                    "combined_features": ef + gf,
                }
            )

        noinput = run_noinput(cell, branch_rows)
        if time_grid is None:
            time_grid = noinput["t"].copy()
        elif not np.allclose(
            time_grid, noinput["t"], rtol=0, atol=1e-12
        ):
            raise RuntimeError("time grid changed across cells")

        cell_receipt = {
            "cell_order": int(model["order"]),
            "species": model["species"],
            "layer": model["layer"],
            "morphology_identifier": model["morphology_identifier"],
            "branches": [],
        }

        for branch in branch_rows:
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
                "electro_features": np.asarray(
                    branch["electro_features"], dtype=float
                ),
                "gross_features": np.asarray(
                    branch["gross_features"], dtype=float
                ),
                "combined_features": np.asarray(
                    branch["combined_features"], dtype=float
                ),
                "G": G.reshape(-1),
                "T": T.reshape(-1),
            }
            rows.append(row)
            cell_receipt["branches"].append(
                {
                    "branch_index": int(branch["branch_index"]),
                    "section": branch["section"],
                    "electro_features": [
                        float(x) for x in branch["electro_features"]
                    ],
                }
            )

        cells.append(cell_receipt)
        print(
            f"[{ci+1:02d}/24] {model['species']:5s} "
            f"{model['morphology_identifier']:>12s} "
            f"electro_cell={cell_ctx['total_dendritic_electrotonic_length']:.3f}"
        )

    if len(rows) != 144:
        raise RuntimeError(f"expected 144 rows, got {len(rows)}")
    return {
        "rows": rows,
        "time": time_grid,
        "cells": cells,
    }


def fit_predict_with_coefficients(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    X_test: np.ndarray,
    Y_test: np.ndarray,
) -> dict:
    Xz, Xtz = standardize(X_train, X_test)
    mean, basis = pca_basis(Y_train, PCA_COMPONENTS)
    scores = (Y_train - mean) @ basis.T

    X1 = np.column_stack([np.ones(len(Xz)), Xz])
    reg = np.eye(X1.shape[1]) * RIDGE_ALPHA
    reg[0, 0] = 0.0
    beta = np.linalg.solve(
        X1.T @ X1 + reg,
        X1.T @ scores,
    )
    pred_scores = np.column_stack(
        [np.ones(len(Xtz)), Xtz]
    ) @ beta
    predicted = mean + pred_scores @ basis

    oracle_scores = (Y_test - mean) @ basis.T
    pca_oracle = mean + oracle_scores @ basis

    mean_pred = np.repeat(mean[None, :], len(X_test), axis=0)

    nearest_idx = []
    for row in Xtz:
        dist = np.sum((Xz - row[None, :]) ** 2, axis=1)
        nearest_idx.append(int(np.argmin(dist)))
    nearest_pred = Y_train[np.asarray(nearest_idx, dtype=int)]

    feature_importance = np.sqrt(
        np.mean(beta[1:, :] ** 2, axis=1)
    )

    return {
        "predicted": predicted,
        "pca_oracle": pca_oracle,
        "mean_pred": mean_pred,
        "nearest_pred": nearest_pred,
        "basis_rank": int(basis.shape[0]),
        "feature_importance": feature_importance,
    }


def loco(rows: list[dict], feature_key: str) -> dict:
    names = {
        "electro_features": ELECTRO_NAMES,
        "gross_features": [f"gross_{i}" for i in range(10)],
        "combined_features": (
            ELECTRO_NAMES + [f"gross_{i}" for i in range(10)]
        ),
    }[feature_key]

    all_scores = {
        "morphology": [],
        "mean": [],
        "nearest": [],
        "pca_oracle": [],
    }
    folds = []
    importance = []

    for fi, heldout in enumerate(
        sorted({int(r["cell_order"]) for r in rows})
    ):
        train = [r for r in rows if int(r["cell_order"]) != heldout]
        test = [r for r in rows if int(r["cell_order"]) == heldout]

        Xtr = row_matrix(train, feature_key)
        Xte = row_matrix(test, feature_key)
        predG = fit_predict_with_coefficients(
            Xtr, row_matrix(train, "G"),
            Xte, row_matrix(test, "G"),
        )
        predT = fit_predict_with_coefficients(
            Xtr, row_matrix(train, "T"),
            Xte, row_matrix(test, "T"),
        )

        importance.append(
            0.5 * (
                predG["feature_importance"]
                + predT["feature_importance"]
            )
        )

        fold = {
            "heldout_cell_order": int(heldout),
            "species": test[0]["species"],
            "morphology_identifier": test[0][
                "morphology_identifier"
            ],
        }
        for model in all_scores:
            scores = score_predictions(
                test, predG, predT, model
            )
            all_scores[model].extend(scores)
            fold[model] = {
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
            f"LOCO {feature_key:18s} [{fi+1:02d}/24] "
            f"{heldout:02d} {test[0]['species']:5s} "
            f"morph={fold['morphology']['median_joint_nrmse']:.3f} "
            f"near={fold['nearest']['median_joint_nrmse']:.3f} "
            f"pca={fold['pca_oracle']['median_joint_nrmse']:.3f}"
        )
        folds.append(fold)

    agg = {
        k: aggregate_scores(v)
        for k, v in all_scores.items()
    }
    wins = sum(
        1 for fold in folds
        if fold["morphology"]["median_joint_nrmse"]
        < fold["nearest"]["median_joint_nrmse"]
    )
    imp = np.mean(np.stack(importance, axis=0), axis=0)
    ranking = sorted(
        [
            {
                "feature": name,
                "mean_rms_standardized_ridge_weight": float(value),
            }
            for name, value in zip(names, imp)
        ],
        key=lambda x: x[
            "mean_rms_standardized_ridge_weight"
        ],
        reverse=True,
    )

    return {
        "feature_key": feature_key,
        "folds": folds,
        "aggregate": agg,
        "morphology_beats_nearest_cells": int(wins),
        "feature_weight_ranking": ranking,
    }


def species_summary(result: dict) -> dict:
    out = {}
    for species in ("rat", "human"):
        vals = [
            f["morphology"]["median_joint_nrmse"]
            for f in result["folds"]
            if f["species"] == species
        ]
        out[species] = {
            "cells": len(vals),
            "median_cell_joint_nrmse": float(np.median(vals)),
        }
    return out


def species_transfer(
    rows: list[dict],
    train_species: str,
    test_species: str,
) -> dict:
    train = [r for r in rows if r["species"] == train_species]
    test = [r for r in rows if r["species"] == test_species]
    Xtr = row_matrix(train, "electro_features")
    Xte = row_matrix(test, "electro_features")

    predG = fit_predict_with_coefficients(
        Xtr, row_matrix(train, "G"),
        Xte, row_matrix(test, "G"),
    )
    predT = fit_predict_with_coefficients(
        Xtr, row_matrix(train, "T"),
        Xte, row_matrix(test, "T"),
    )
    scores = score_predictions(
        test, predG, predT, "morphology"
    )
    return aggregate_scores(scores)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/electrotonic_chart.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source not pinned")

    setup_neuron(fci_root)
    measured = collect_panel(
        fci_root, load_panel(fci_root)
    )
    rows = measured["rows"]

    electro = loco(rows, "electro_features")
    gross = loco(rows, "gross_features")
    combined = loco(rows, "combined_features")

    e = electro["aggregate"]["morphology"]
    near = electro["aggregate"]["nearest"]
    pca = electro["aggregate"]["pca_oracle"]
    previous_gross = 0.3522

    ratio_previous = (
        e["median_joint_nrmse"] / previous_gross
    )
    ratio_nearest = (
        e["median_joint_nrmse"]
        / (near["median_joint_nrmse"] + 1e-30)
    )

    strong = e["median_joint_nrmse"] <= 0.20
    improve = (
        e["median_joint_nrmse"] <= 0.30
        and e["median_G_nrmse"] <= 0.30
        and e["median_T_nrmse"] <= 0.32
        and ratio_previous <= 0.85
        and ratio_nearest <= 0.90
        and electro["morphology_beats_nearest_cells"] >= 16
        and pca["median_joint_nrmse"] <= 0.10
    )

    if strong and improve:
        classification = (
            "CROSS_CELL_OPERATOR_PREDICTABLE_FROM_ELECTROTONIC_MORPHOLOGY"
        )
    elif improve:
        classification = (
            "ELECTROTONIC_CHART_IMPROVES_CROSS_CELL_OPERATOR_MAP"
        )
    else:
        classification = (
            "OPERATOR_LOW_DIMENSIONAL_ELECTROTONIC_CHART_STILL_INSUFFICIENT"
        )

    summary = {
        "object": (
            "cross-cell operator prediction from morphology-derived "
            "electrotonic and Rall-style coordinates"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "cells": 24,
            "branches_per_cell": N_BRANCHES,
            "features": ELECTRO_NAMES,
            "matched_passive": {
                "Ra_ohm_cm": RA_OHM_CM,
                "Rm_ohm_cm2": RM_OHM_CM2,
            },
            "pca_components": PCA_COMPONENTS,
            "ridge_alpha": RIDGE_ALPHA,
            "heldout_unit": "whole cell",
            "previous_gross_map_median_joint_nrmse": previous_gross,
            "thresholds_locked_before_run": {
                "joint_nrmse_max": 0.30,
                "G_nrmse_max": 0.30,
                "T_nrmse_max": 0.32,
                "to_previous_gross_ratio_max": 0.85,
                "to_nearest_ratio_max": 0.90,
                "beats_nearest_cells_min": 16,
                "pca_oracle_max": 0.10,
                "strong_joint_nrmse_max": 0.20,
            },
        },
        "electrotonic": electro,
        "gross_rerun": gross,
        "combined": combined,
        "primary_ratios": {
            "electrotonic_to_previous_gross": float(
                ratio_previous
            ),
            "electrotonic_to_nearest": float(
                ratio_nearest
            ),
        },
        "species": species_summary(electro),
        "species_transfer": {
            "rat_to_human": species_transfer(
                rows, "rat", "human"
            ),
            "human_to_rat": species_transfer(
                rows, "human", "rat"
            ),
        },
        "classification": classification,
        "stopping_line": (
            "Do not add neural nets, polynomial features, impedance "
            "measurements, species labels, target-cell fitting or delete "
            "outliers after this result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori electrotonic-chart audit")
    print()
    print(
        "electrotonic joint NRMSE:       "
        f"{e['median_joint_nrmse']:.4f}"
    )
    print(
        "  G NRMSE:                      "
        f"{e['median_G_nrmse']:.4f}"
    )
    print(
        "  T NRMSE:                      "
        f"{e['median_T_nrmse']:.4f}"
    )
    print(
        "electro nearest NRMSE:           "
        f"{near['median_joint_nrmse']:.4f}"
    )
    print(
        "PCA oracle NRMSE:                "
        f"{pca['median_joint_nrmse']:.4f}"
    )
    print(
        "electro / previous gross:        "
        f"{ratio_previous:.4f}"
    )
    print(
        "electro / nearest:               "
        f"{ratio_nearest:.4f}"
    )
    print(
        "beats nearest cells:             "
        f"{electro['morphology_beats_nearest_cells']} / 24"
    )
    print(
        "gross rerun joint NRMSE:         "
        f"{gross['aggregate']['morphology']['median_joint_nrmse']:.4f}"
    )
    print(
        "combined joint NRMSE:            "
        f"{combined['aggregate']['morphology']['median_joint_nrmse']:.4f}"
    )
    print(
        "rat median cell joint:           "
        f"{species_summary(electro)['rat']['median_cell_joint_nrmse']:.4f}"
    )
    print(
        "human median cell joint:         "
        f"{species_summary(electro)['human']['median_cell_joint_nrmse']:.4f}"
    )
    print(
        "rat -> human joint:              "
        f"{summary['species_transfer']['rat_to_human']['median_joint_nrmse']:.4f}"
    )
    print(
        "human -> rat joint:              "
        f"{summary['species_transfer']['human_to_rat']['median_joint_nrmse']:.4f}"
    )
    print("top electrotonic features:")
    for row in electro["feature_weight_ranking"][:6]:
        print(
            f"  {row['feature']:44s} "
            f"{row['mean_rms_standardized_ridge_weight']:.4f}"
        )
    print(f"classification: {classification}")

    assert len(rows) == 144


if __name__ == "__main__":
    main()
