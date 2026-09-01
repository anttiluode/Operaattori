from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from gate17_superposition_attack import (
    FCI_COMMIT,
    MODEL_REL,
    V_INIT_MV,
    compact_midpoint_rows,
    configure_active,
    dendritic_rows,
    git_head,
    import_builder,
    passive_site_features,
    positive_auc,
    run_trace,
    section_groups,
    settle_baselines,
)


DOSES = (4.0, 8.0, 12.0)
SPLITS = ((4.0, 12.0), (8.0, 8.0), (12.0, 4.0))


def run_weighted_pair(
    cell,
    syn_df,
    sites_a: np.ndarray,
    sites_b: np.ndarray,
    multiplicity_a: float,
    multiplicity_b: float,
    condition: str,
    baseline_by_row: dict[int, float],
    *,
    event_ms: float = 60.0,
    tstop_ms: float = 160.0,
) -> dict:
    from neuron import h

    sites_a = np.asarray(sites_a, dtype=int)
    sites_b = np.asarray(sites_b, dtype=int)
    record_rows = np.concatenate([sites_a, sites_b])

    configure_active(
        syn_df,
        sites_a,
        float(multiplicity_a),
        condition,
        baseline_by_row,
    )
    configure_active(
        syn_df,
        sites_b,
        float(multiplicity_b),
        condition,
        baseline_by_row,
    )

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)

    h.dt = 0.025
    h.finitialize(V_INIT_MV)
    h.fcurrent()
    for i in record_rows:
        syn_df.iloc[int(i)]["exc_netcons"].event(float(event_ms))
    h.continuerun(float(tstop_ms))

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    if not np.all(np.isfinite(soma)):
        raise FloatingPointError("non-finite Gate-21 soma trace")

    pre = (t >= event_ms - 10.0) & (t < event_ms - 1.0)
    post = (t >= event_ms) & (t <= event_ms + 90.0)
    soma_base = float(np.median(soma[pre]))
    depol = soma[post] - soma_base

    absolute = depol + soma_base
    above = absolute >= -20.0
    spike_crossings = (
        int(np.sum((~above[:-1]) & above[1:])) if len(above) > 1 else 0
    )
    large_depol = bool(np.max(depol) >= 40.0)

    return {
        "t": t[post],
        "soma_baseline_mV": soma_base,
        "soma_depol": depol,
        "conservative_soma_spike_crossings": spike_crossings,
        "large_soma_depolarization_guard": large_depol,
    }


def nrmse(actual: np.ndarray, predicted: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    denom = float(np.sqrt(np.mean(actual * actual))) + 1e-30
    return float(np.sqrt(np.mean((predicted - actual) ** 2)) / denom)


def log_signature(values: list[float]) -> np.ndarray:
    arr = np.maximum(np.asarray(values, dtype=float), 1e-30)
    logs = np.log(arr)
    return logs - float(np.mean(logs))


def sig_rmse(actual: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.sqrt(np.mean((actual - predicted) ** 2)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=6)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/gate21/gate21_branch_subunit_computation.json"
        ),
    )
    args = ap.parse_args()

    if args.branches != 6 or args.sites != 3:
        raise ValueError("Gate 21 is locked to six branches x three sites")
    if git_head(args.fci_root.resolve()) != FCI_COMMIT:
        raise RuntimeError("FCI checkout is not pinned")

    builder = import_builder(args.fci_root.resolve())
    cell, syn_df = builder.create_cell(
        path=str(args.fci_root.resolve() / MODEL_REL) + "/"
    )

    rows = dendritic_rows(syn_df)
    groups = section_groups(syn_df, rows)
    _zinput, _transfer, _run_ids, _name_to_id = passive_site_features(
        cell, syn_df, rows
    )
    baseline_by_row = settle_baselines(syn_df, rows)

    candidates = []
    for name, row_ids in groups.items():
        clustered, span = compact_midpoint_rows(
            syn_df,
            row_ids,
            args.sites,
            args.cluster_span_um,
        )
        if len(clustered) != args.sites:
            continue
        clustered = np.asarray(
            sorted(
                clustered.tolist(),
                key=lambda i: float(syn_df.iloc[int(i)]["segments"].x),
            ),
            dtype=int,
        )
        sec = syn_df.iloc[int(clustered[0])]["segments"].sec
        candidates.append(
            {
                "section": name,
                "section_length_um": float(sec.L),
                "sites": clustered,
                "span_um": float(span),
            }
        )

    candidates.sort(
        reverse=True,
        key=lambda x: x["section_length_um"],
    )
    candidates = candidates[:6]
    if len(candidates) != 6:
        raise RuntimeError("Gate 21 requires the six Gate-16/20 sections")

    conditions = ("human", "gamma062_restmatched")

    # Precompute the two nested physical rulers:
    #   site_sum[b,d]   = sum of three independent single-site traces
    #   branch[b,d]     = the complete nonlinear three-site branch-alone trace
    basis: dict[str, dict[int, dict[float, dict]]] = {}

    for condition in conditions:
        basis[condition] = {}
        for bi, branch in enumerate(candidates):
            basis[condition][bi] = {}
            for dose in DOSES:
                branch_trace = run_trace(
                    cell,
                    syn_df,
                    branch["sites"],
                    branch["sites"],
                    dose,
                    condition,
                    baseline_by_row,
                )
                singles = [
                    run_trace(
                        cell,
                        syn_df,
                        branch["sites"],
                        np.asarray([site], dtype=int),
                        dose,
                        condition,
                        baseline_by_row,
                    )
                    for site in branch["sites"]
                ]
                t = branch_trace["t"]
                if any(
                    not np.allclose(x["t"], t, rtol=0, atol=1e-12)
                    for x in singles
                ):
                    raise RuntimeError("Gate-21 basis time grids differ")
                site_sum = np.sum(
                    np.stack([x["soma_depol"] for x in singles], axis=0),
                    axis=0,
                )
                basis[condition][bi][dose] = {
                    "t": t,
                    "branch_trace": branch_trace["soma_depol"],
                    "site_sum_trace": site_sum,
                    "branch_auc": positive_auc(
                        branch_trace["soma_depol"], t
                    ),
                    "site_sum_auc": positive_auc(site_sum, t),
                }
                print(
                    f"basis {condition:22s} branch={bi} dose={dose:4.0f} "
                    f"branchAUC={basis[condition][bi][dose]['branch_auc']:.4f} "
                    f"siteAUC={basis[condition][bi][dose]['site_sum_auc']:.4f}"
                )

    pair_rows = []
    for pi, (ia, ib) in enumerate(
        itertools.combinations(range(len(candidates)), 2)
    ):
        prow = {
            "pair_index": int(pi),
            "branch_a": int(ia),
            "branch_b": int(ib),
            "section_a": candidates[ia]["section"],
            "section_b": candidates[ib]["section"],
            "conditions": {},
        }

        for condition in conditions:
            split_rows = []
            actual_aucs: list[float] = []
            site_aucs: list[float] = []
            subunit_aucs: list[float] = []
            site_trace_errors: list[float] = []
            subunit_trace_errors: list[float] = []
            any_spike = False

            for ma, mb in SPLITS:
                actual = run_weighted_pair(
                    cell,
                    syn_df,
                    candidates[ia]["sites"],
                    candidates[ib]["sites"],
                    ma,
                    mb,
                    condition,
                    baseline_by_row,
                )
                ta = basis[condition][ia][ma]["t"]
                tb = basis[condition][ib][mb]["t"]
                if (
                    not np.allclose(actual["t"], ta, rtol=0, atol=1e-12)
                    or not np.allclose(actual["t"], tb, rtol=0, atol=1e-12)
                ):
                    raise RuntimeError("Gate-21 pair/basis time grids differ")

                site_pred = (
                    basis[condition][ia][ma]["site_sum_trace"]
                    + basis[condition][ib][mb]["site_sum_trace"]
                )
                subunit_pred = (
                    basis[condition][ia][ma]["branch_trace"]
                    + basis[condition][ib][mb]["branch_trace"]
                )

                actual_auc = positive_auc(
                    actual["soma_depol"], actual["t"]
                )
                site_auc = positive_auc(site_pred, actual["t"])
                subunit_auc = positive_auc(subunit_pred, actual["t"])

                serr = nrmse(actual["soma_depol"], site_pred)
                berr = nrmse(actual["soma_depol"], subunit_pred)
                spike = bool(
                    actual["conservative_soma_spike_crossings"] > 0
                    or actual["large_soma_depolarization_guard"]
                )
                any_spike = any_spike or spike

                actual_aucs.append(actual_auc)
                site_aucs.append(site_auc)
                subunit_aucs.append(subunit_auc)
                site_trace_errors.append(serr)
                subunit_trace_errors.append(berr)

                split_rows.append(
                    {
                        "multiplicity_a_per_site": ma,
                        "multiplicity_b_per_site": mb,
                        "virtual_synapses_a": int(round(3.0 * ma)),
                        "virtual_synapses_b": int(round(3.0 * mb)),
                        "total_virtual_synapses": int(
                            round(3.0 * (ma + mb))
                        ),
                        "actual_soma_auc_mV_ms": actual_auc,
                        "independent_site_pred_auc_mV_ms": site_auc,
                        "nonlinear_subunit_pred_auc_mV_ms": subunit_auc,
                        "independent_site_trace_nrmse": serr,
                        "nonlinear_subunit_trace_nrmse": berr,
                        "spike_guard": spike,
                    }
                )

            actual_sig = log_signature(actual_aucs)
            site_sig = log_signature(site_aucs)
            subunit_sig = log_signature(subunit_aucs)
            zero_sig = np.zeros_like(actual_sig)

            total_err = sig_rmse(actual_sig, zero_sig)
            site_err = sig_rmse(actual_sig, site_sig)
            subunit_err = sig_rmse(actual_sig, subunit_sig)
            range_factor = float(
                max(actual_aucs) / (min(actual_aucs) + 1e-30)
            )

            prow["conditions"][condition] = {
                "splits": split_rows,
                "actual_centered_log_auc_signature": actual_sig.tolist(),
                "independent_site_signature": site_sig.tolist(),
                "nonlinear_subunit_signature": subunit_sig.tolist(),
                "equal_budget_auc_range_factor": range_factor,
                "total_input_signature_rmse": total_err,
                "independent_site_signature_rmse": site_err,
                "nonlinear_subunit_signature_rmse": subunit_err,
                "median_independent_site_trace_nrmse": float(
                    np.median(site_trace_errors)
                ),
                "median_nonlinear_subunit_trace_nrmse": float(
                    np.median(subunit_trace_errors)
                ),
                "nonlinear_subunit_beats_site_signature": bool(
                    subunit_err < site_err
                ),
                "spike_guard": bool(any_spike),
            }

        pair_rows.append(prow)
        hrow = prow["conditions"]["human"]
        print(
            f"pair [{pi+1:02d}/15] {ia}-{ib} "
            f"range={hrow['equal_budget_auc_range_factor']:.4f}x "
            f"sig total/site/subunit="
            f"{hrow['total_input_signature_rmse']:.4f}/"
            f"{hrow['independent_site_signature_rmse']:.4f}/"
            f"{hrow['nonlinear_subunit_signature_rmse']:.4f}"
        )

    aggregate = {}
    for condition in conditions:
        rows_c = [p["conditions"][condition] for p in pair_rows]
        ranges = np.asarray(
            [x["equal_budget_auc_range_factor"] for x in rows_c],
            dtype=float,
        )
        total_err = np.asarray(
            [x["total_input_signature_rmse"] for x in rows_c],
            dtype=float,
        )
        site_err = np.asarray(
            [x["independent_site_signature_rmse"] for x in rows_c],
            dtype=float,
        )
        subunit_err = np.asarray(
            [x["nonlinear_subunit_signature_rmse"] for x in rows_c],
            dtype=float,
        )
        site_trace = np.asarray(
            [x["median_independent_site_trace_nrmse"] for x in rows_c],
            dtype=float,
        )
        sub_trace = np.asarray(
            [x["median_nonlinear_subunit_trace_nrmse"] for x in rows_c],
            dtype=float,
        )
        spikes = np.asarray(
            [x["spike_guard"] for x in rows_c],
            dtype=bool,
        )

        aggregate[condition] = {
            "median_equal_budget_auc_range_factor": float(np.median(ranges)),
            "fraction_pairs_range_at_least_1p05": float(
                np.mean(ranges >= 1.05)
            ),
            "median_total_input_signature_rmse": float(
                np.median(total_err)
            ),
            "median_independent_site_signature_rmse": float(
                np.median(site_err)
            ),
            "median_nonlinear_subunit_signature_rmse": float(
                np.median(subunit_err)
            ),
            "site_to_subunit_signature_error_ratio": float(
                np.median(site_err) / (np.median(subunit_err) + 1e-30)
            ),
            "fraction_pairs_subunit_beats_site": float(
                np.mean(subunit_err < site_err)
            ),
            "fraction_pairs_subunit_beats_total": float(
                np.mean(subunit_err < total_err)
            ),
            "median_independent_site_trace_nrmse": float(
                np.median(site_trace)
            ),
            "median_nonlinear_subunit_trace_nrmse": float(
                np.median(sub_trace)
            ),
            "fraction_pairs_spike_guard": float(np.mean(spikes)),
        }

    h = aggregate["human"]
    if h["fraction_pairs_spike_guard"] > 0.0:
        classification = "EQUAL_BUDGET_ASSAY_SPIKING_CONFOUNDED"
        interpretation = (
            "At least one equal-budget redistribution crossed the somatic "
            "spike guard, so the branch-factorization comparison is not a "
            "clean subthreshold assay."
        )
    elif (
        h["median_equal_budget_auc_range_factor"] < 1.05
        or h["fraction_pairs_range_at_least_1p05"] < 0.50
    ):
        classification = "EQUAL_BUDGET_REDISTRIBUTION_WEAK"
        interpretation = (
            "Redistributing the fixed 48-synapse budget across the same two "
            "branch compartments does not robustly create a >=5% somatic "
            "output signature."
        )
    elif (
        h["site_to_subunit_signature_error_ratio"] >= 2.0
        and h["fraction_pairs_subunit_beats_site"] >= 0.80
        and h["fraction_pairs_subunit_beats_total"] >= 0.80
        and h["median_nonlinear_subunit_trace_nrmse"] <= 0.10
    ):
        classification = (
            "NONLINEAR_BRANCH_SUBUNITS_PREDICT_EQUAL_BUDGET_OUTPUTS"
        )
        interpretation = (
            "Across fixed branch pairs and fixed total synapse count, the "
            "sum of the two measured nonlinear branch-alone responses predicts "
            "the redistribution-dependent somatic signature substantially "
            "better than the sum of independent sites. The scaffold therefore "
            "admits a useful nonlinear-subunit factorization at these doses."
        )
    elif (
        h["median_independent_site_signature_rmse"]
        <= 1.10 * h["median_nonlinear_subunit_signature_rmse"]
    ):
        classification = "INDEPENDENT_SITE_TRANSFER_SUFFICIENT"
        interpretation = (
            "The nonlinear branch basis does not materially improve the "
            "equal-budget redistribution signature over independent site "
            "transfer. The Gate-20 compartment label does not buy a simpler "
            "somatic computation model here."
        )
    else:
        classification = "NO_CLEAN_BRANCH_SUBUNIT_FACTORISATION"
        interpretation = (
            "Equal-budget redistribution changes the somatic output, but the "
            "preregistered two-branch nonlinear-subunit sum does not explain "
            "the signatures cleanly enough to earn the factorized model."
        )

    summary = {
        "gate": 21,
        "object": (
            "equal-total-input redistribution across measured nonlinear "
            "branch subunits"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "branches": 6,
            "sites_per_branch": 3,
            "all_branch_pairs": 15,
            "per_site_doses": list(DOSES),
            "equal_budget_splits": [list(x) for x in SPLITS],
            "total_virtual_synapses_each_pattern": 48,
            "same_two_branches_and_same_six_sites_within_each_pair": True,
            "primary_output": "somatic positive depolarization AUC",
            "primary_signature": (
                "centered log AUC across 4+12, 8+8, 12+4 per-site "
                "multiplicity splits"
            ),
            "total_input_ruler": "zero centered redistribution signature",
            "independent_site_ruler": (
                "sum of six matching single-site soma traces"
            ),
            "nonlinear_subunit_ruler": (
                "sum of the two complete branch-alone soma traces at the "
                "exact two doses; no pair response is used"
            ),
            "positive_thresholds_locked_before_run": {
                "median_range_factor_min": 1.05,
                "fraction_pairs_range_min": 0.50,
                "site_to_subunit_signature_error_ratio_min": 2.0,
                "fraction_pairs_subunit_beats_site_min": 0.80,
                "fraction_pairs_subunit_beats_total_min": 0.80,
                "median_subunit_trace_nrmse_max": 0.10,
            },
            "conditions": list(conditions),
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "pairs": pair_rows,
        "stopping_line": (
            "A positive Gate 21 earns only a held-out-dose/generalization "
            "attack on the branch transfer law. Do not add growth yet. A "
            "negative result is not to be rescued by scanning dose splits."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori Gate 21 — equal-budget branch-subunit computation")
    print()
    for condition in conditions:
        a = aggregate[condition]
        print(condition.upper())
        print(
            "  median equal-budget range:          "
            f"{a['median_equal_budget_auc_range_factor']:.4f}x"
        )
        print(
            "  pairs range >=1.05:                "
            f"{a['fraction_pairs_range_at_least_1p05']:.3f}"
        )
        print(
            "  signature RMSE total/site/subunit: "
            f"{a['median_total_input_signature_rmse']:.4f} / "
            f"{a['median_independent_site_signature_rmse']:.4f} / "
            f"{a['median_nonlinear_subunit_signature_rmse']:.4f}"
        )
        print(
            "  site/subunit error ratio:          "
            f"{a['site_to_subunit_signature_error_ratio']:.3f}x"
        )
        print(
            "  subunit beats site / total:        "
            f"{a['fraction_pairs_subunit_beats_site']:.3f} / "
            f"{a['fraction_pairs_subunit_beats_total']:.3f}"
        )
        print(
            "  trace NRMSE site/subunit:           "
            f"{a['median_independent_site_trace_nrmse']:.4f} / "
            f"{a['median_nonlinear_subunit_trace_nrmse']:.4f}"
        )
        print(
            "  spike guard:                       "
            f"{a['fraction_pairs_spike_guard']:.3f}"
        )
        print()

    print(f"classification: {classification}")
    print(interpretation)

    assert len(pair_rows) == 15
    assert all(
        np.all(np.isfinite([
            p["conditions"][c]["equal_budget_auc_range_factor"],
            p["conditions"][c]["total_input_signature_rmse"],
            p["conditions"][c]["independent_site_signature_rmse"],
            p["conditions"][c]["nonlinear_subunit_signature_rmse"],
        ]))
        for p in pair_rows
        for c in conditions
    )
    assert classification in {
        "EQUAL_BUDGET_ASSAY_SPIKING_CONFOUNDED",
        "EQUAL_BUDGET_REDISTRIBUTION_WEAK",
        "NONLINEAR_BRANCH_SUBUNITS_PREDICT_EQUAL_BUDGET_OUTPUTS",
        "INDEPENDENT_SITE_TRANSFER_SUFFICIENT",
        "NO_CLEAN_BRANCH_SUBUNIT_FACTORISATION",
    }


if __name__ == "__main__":
    main()
