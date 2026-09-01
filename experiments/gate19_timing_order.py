from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from gate17_superposition_attack import (
    FCI_COMMIT,
    MODEL_REL,
    compact_midpoint_rows,
    configure_active,
    dendritic_rows,
    git_head,
    import_builder,
    passive_site_features,
    positive_auc,
    section_groups,
    settle_baselines,
)
from operaattori.compartment_match import greedy_dispersed_match


def run_timed_trace(
    cell,
    syn_df,
    record_rows: np.ndarray,
    active_rows: np.ndarray,
    event_times_ms: np.ndarray,
    multiplicity: float,
    condition: str,
    baseline_by_row: dict[int, float],
    *,
    tstop_ms: float,
    analysis_first_ms: float | None = None,
    analysis_last_ms: float | None = None,
) -> dict:
    from neuron import h

    record_rows = np.asarray(record_rows, dtype=int)
    active_rows = np.asarray(active_rows, dtype=int)
    event_times_ms = np.asarray(event_times_ms, dtype=float)

    if len(active_rows) != len(event_times_ms):
        raise ValueError("one event time is required for each active site")

    configure_active(
        syn_df,
        active_rows,
        float(multiplicity),
        condition,
        baseline_by_row,
    )

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
    local_vecs = [
        h.Vector().record(syn_df.iloc[int(i)]["segments"]._ref_v)
        for i in record_rows
    ]

    h.dt = 0.025
    h.finitialize(-76.0)
    h.fcurrent()

    for row_id, event_ms in zip(active_rows, event_times_ms):
        syn_df.iloc[int(row_id)]["exc_netcons"].event(float(event_ms))

    h.continuerun(float(tstop_ms))

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    local = np.stack(
        [np.asarray(v, dtype=float) for v in local_vecs],
        axis=0,
    )

    if not np.all(np.isfinite(soma)) or not np.all(np.isfinite(local)):
        raise FloatingPointError("non-finite Gate-19 trace")

    first_event = (
        float(np.min(event_times_ms))
        if analysis_first_ms is None
        else float(analysis_first_ms)
    )
    last_event = (
        float(np.max(event_times_ms))
        if analysis_last_ms is None
        else float(analysis_last_ms)
    )
    pre = (t >= first_event - 10.0) & (t < first_event - 1.0)
    post = (t >= first_event) & (t <= last_event + 90.0)

    if not np.any(pre) or not np.any(post):
        raise RuntimeError("Gate-19 baseline/post windows missing")

    soma_base = float(np.median(soma[pre]))
    local_base = np.median(local[:, pre], axis=1)

    return {
        "t": t[post],
        "soma_depol": soma[post] - soma_base,
        "local_depol": local[:, post] - local_base[:, None],
    }


def timed_interaction_metrics(
    cell,
    syn_df,
    sites: np.ndarray,
    event_times_ms: np.ndarray,
    multiplicity: float,
    condition: str,
    baseline_by_row: dict[int, float],
) -> dict:
    sites = np.asarray(sites, dtype=int)
    event_times_ms = np.asarray(event_times_ms, dtype=float)
    analysis_first = float(np.min(event_times_ms))
    analysis_last = float(np.max(event_times_ms))
    tstop = float(analysis_last + 100.0)

    together = run_timed_trace(
        cell,
        syn_df,
        sites,
        sites,
        event_times_ms,
        multiplicity,
        condition,
        baseline_by_row,
        tstop_ms=tstop,
        analysis_first_ms=analysis_first,
        analysis_last_ms=analysis_last,
    )

    singles = []
    for site, event_ms in zip(sites, event_times_ms):
        singles.append(
            run_timed_trace(
                cell,
                syn_df,
                sites,
                np.asarray([site], dtype=int),
                np.asarray([event_ms], dtype=float),
                multiplicity,
                condition,
                baseline_by_row,
                tstop_ms=tstop,
                analysis_first_ms=analysis_first,
                analysis_last_ms=analysis_last,
            )
        )

    t = together["t"]
    for single in singles:
        if not np.allclose(single["t"], t, rtol=0, atol=1e-12):
            raise RuntimeError("Gate-19 single-site time grids differ")

    predicted_local = np.sum(
        np.stack([x["local_depol"] for x in singles], axis=0),
        axis=0,
    )
    predicted_soma = np.sum(
        np.stack([x["soma_depol"] for x in singles], axis=0),
        axis=0,
    )

    actual_local_mean = np.mean(together["local_depol"], axis=0)
    predicted_local_mean = np.mean(predicted_local, axis=0)

    actual_local_auc = positive_auc(actual_local_mean, t)
    predicted_local_auc = positive_auc(predicted_local_mean, t)
    actual_soma_auc = positive_auc(together["soma_depol"], t)
    predicted_soma_auc = positive_auc(predicted_soma, t)

    return {
        "actual_local_mean_auc_mV_ms": float(actual_local_auc),
        "independent_sum_local_mean_auc_mV_ms": float(predicted_local_auc),
        "local_interaction_ratio": float(
            actual_local_auc / (predicted_local_auc + 1e-30)
        ),
        "actual_soma_auc_mV_ms": float(actual_soma_auc),
        "independent_sum_soma_auc_mV_ms": float(predicted_soma_auc),
        "soma_interaction_ratio": float(
            actual_soma_auc / (predicted_soma_auc + 1e-30)
        ),
        "event_times_ms": event_times_ms.tolist(),
    }


def order_metrics(forward: dict, reverse: dict) -> dict:
    fi = float(forward["local_interaction_ratio"])
    ri = float(reverse["local_interaction_ratio"])
    signed_log = float(math.log((fi + 1e-30) / (ri + 1e-30)))
    return {
        "forward_interaction_ratio": fi,
        "reverse_interaction_ratio": ri,
        "signed_log_forward_over_reverse": signed_log,
        "order_magnitude_log": abs(signed_log),
        "order_magnitude_fraction": float(math.exp(abs(signed_log)) - 1.0),
        "preferred_direction": (
            "proximal_to_distal" if signed_log > 0
            else "distal_to_proximal" if signed_log < 0
            else "tie"
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=6)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument("--multiplicity", type=float, default=16.0)
    ap.add_argument("--lag-ms", type=float, default=4.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/gate19/gate19_timing_order.json"),
    )
    args = ap.parse_args()

    if args.sites != 3:
        raise ValueError("Gate 19 protocol is locked to three sites")
    if abs(float(args.lag_ms) - 4.0) > 1e-12:
        raise ValueError("Gate 19 primary inter-site lag is locked to 4 ms")
    if git_head(args.fci_root.resolve()) != FCI_COMMIT:
        raise RuntimeError("FCI checkout is not pinned")

    builder = import_builder(args.fci_root.resolve())
    cell, syn_df = builder.create_cell(
        path=str(args.fci_root.resolve() / MODEL_REL) + "/"
    )

    rows = dendritic_rows(syn_df)
    groups = section_groups(syn_df, rows)
    zinput, soma_transfer, run_ids, name_to_id = passive_site_features(
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
        candidates.append((float(sec.L), name, clustered, span))

    candidates.sort(reverse=True, key=lambda x: x[0])
    candidates = candidates[: min(args.branches, len(candidates))]
    if len(candidates) < 4:
        raise RuntimeError("too few compact sections for Gate 19")

    t0 = 60.0
    forward_times = np.asarray(
        [t0, t0 + args.lag_ms, t0 + 2.0 * args.lag_ms],
        dtype=float,
    )
    reverse_times = forward_times[::-1].copy()

    threshold_log = float(math.log(1.05))
    result_rows = []

    for bi, (section_length, name, clustered, span) in enumerate(candidates):
        rid = int(name_to_id[name])
        dispersed, match = greedy_dispersed_match(
            clustered,
            rid,
            rows,
            run_ids,
            zinput,
            soma_transfer,
            min_distinct_runs=3,
        )
        dispersed = np.asarray(dispersed, dtype=int)

        conditions = {}
        for arrangement, sites in (
            ("clustered", clustered),
            ("dispersed", dispersed),
        ):
            conditions[arrangement] = {}
            for condition in ("human", "gamma062_restmatched"):
                forward = timed_interaction_metrics(
                    cell,
                    syn_df,
                    sites,
                    forward_times,
                    float(args.multiplicity),
                    condition,
                    baseline_by_row,
                )
                reverse = timed_interaction_metrics(
                    cell,
                    syn_df,
                    sites,
                    reverse_times,
                    float(args.multiplicity),
                    condition,
                    baseline_by_row,
                )
                conditions[arrangement][condition] = {
                    "forward": forward,
                    "reverse": reverse,
                    "order": order_metrics(forward, reverse),
                }

        ch = conditions["clustered"]["human"]["order"]
        dh = conditions["dispersed"]["human"]["order"]
        cr = conditions["clustered"]["gamma062_restmatched"]["order"]
        dr = conditions["dispersed"]["gamma062_restmatched"]["order"]

        human_excess = float(
            ch["order_magnitude_log"] - dh["order_magnitude_log"]
        )
        rest_excess = float(
            cr["order_magnitude_log"] - dr["order_magnitude_log"]
        )
        gamma_specific_excess = float(human_excess - rest_excess)

        result = {
            "branch_index": int(bi),
            "section": name,
            "section_length_um": float(section_length),
            "cluster_span_um": float(span),
            "clustered_rows_proximal_to_distal": clustered.tolist(),
            "clustered_x": [
                float(syn_df.iloc[int(i)]["segments"].x)
                for i in clustered
            ],
            "dispersed_rows_mapped_to_cluster_order": dispersed.tolist(),
            "passive_match": match,
            "conditions": conditions,
            "human_compact_minus_dispersed_order_magnitude_log": human_excess,
            "restmatched_compact_minus_dispersed_order_magnitude_log": rest_excess,
            "gamma_specific_order_locality_excess_log": gamma_specific_excess,
            "human_compact_order_above_5pct": bool(
                ch["order_magnitude_log"] >= threshold_log
            ),
            "human_compact_more_order_sensitive_than_dispersed": bool(
                human_excess > 0
            ),
        }
        result_rows.append(result)

        print(
            f"[{bi+1:02d}/{len(candidates):02d}] {name} "
            f"span={span:.1f}um "
            f"Hcompact={ch['order_magnitude_fraction']:.3f} "
            f"Hspread={dh['order_magnitude_fraction']:.3f} "
            f"Rcompact={cr['order_magnitude_fraction']:.3f} "
            f"Rspread={dr['order_magnitude_fraction']:.3f} "
            f"H-excess-log={human_excess:.4f} "
            f"gamma-excess-log={gamma_specific_excess:.4f}"
        )

    human_compact_mag = np.asarray(
        [
            r["conditions"]["clustered"]["human"]["order"][
                "order_magnitude_log"
            ]
            for r in result_rows
        ],
        dtype=float,
    )
    human_spread_mag = np.asarray(
        [
            r["conditions"]["dispersed"]["human"]["order"][
                "order_magnitude_log"
            ]
            for r in result_rows
        ],
        dtype=float,
    )
    rest_compact_mag = np.asarray(
        [
            r["conditions"]["clustered"]["gamma062_restmatched"]["order"][
                "order_magnitude_log"
            ]
            for r in result_rows
        ],
        dtype=float,
    )
    rest_spread_mag = np.asarray(
        [
            r["conditions"]["dispersed"]["gamma062_restmatched"]["order"][
                "order_magnitude_log"
            ]
            for r in result_rows
        ],
        dtype=float,
    )

    human_excess = np.asarray(
        [
            r["human_compact_minus_dispersed_order_magnitude_log"]
            for r in result_rows
        ],
        dtype=float,
    )
    rest_excess = np.asarray(
        [
            r["restmatched_compact_minus_dispersed_order_magnitude_log"]
            for r in result_rows
        ],
        dtype=float,
    )
    gamma_excess = np.asarray(
        [
            r["gamma_specific_order_locality_excess_log"]
            for r in result_rows
        ],
        dtype=float,
    )
    match_r = np.asarray(
        [r["passive_match"]["median_z_ratio_factor"] for r in result_rows],
        dtype=float,
    )
    match_t = np.asarray(
        [
            r["passive_match"]["median_transfer_ratio_factor"]
            for r in result_rows
        ],
        dtype=float,
    )

    human_preference = [
        r["conditions"]["clustered"]["human"]["order"]["preferred_direction"]
        for r in result_rows
    ]

    aggregate = {
        "branches": int(len(result_rows)),
        "sites": int(args.sites),
        "multiplicity_per_site": float(args.multiplicity),
        "total_virtual_synapses": int(
            round(args.multiplicity * args.sites)
        ),
        "lag_ms": float(args.lag_ms),
        "event_times_forward_ms": forward_times.tolist(),
        "event_times_reverse_ms": reverse_times.tolist(),
        "median_passive_Rinput_match_factor": float(np.median(match_r)),
        "median_passive_soma_transfer_match_factor": float(np.median(match_t)),
        "human": {
            "median_compact_order_magnitude_fraction": float(
                np.exp(np.median(human_compact_mag)) - 1.0
            ),
            "median_dispersed_order_magnitude_fraction": float(
                np.exp(np.median(human_spread_mag)) - 1.0
            ),
            "fraction_compact_order_above_5pct": float(
                np.mean(human_compact_mag >= threshold_log)
            ),
            "median_compact_minus_dispersed_order_magnitude_log": float(
                np.median(human_excess)
            ),
            "fraction_compact_more_order_sensitive_than_dispersed": float(
                np.mean(human_excess > 0)
            ),
            "proximal_to_distal_preferred_fraction": float(
                np.mean(
                    np.asarray(human_preference, dtype=object)
                    == "proximal_to_distal"
                )
            ),
        },
        "restmatched_gamma062": {
            "median_compact_order_magnitude_fraction": float(
                np.exp(np.median(rest_compact_mag)) - 1.0
            ),
            "median_dispersed_order_magnitude_fraction": float(
                np.exp(np.median(rest_spread_mag)) - 1.0
            ),
            "median_compact_minus_dispersed_order_magnitude_log": float(
                np.median(rest_excess)
            ),
        },
        "gamma_specific": {
            "median_order_locality_excess_log": float(
                np.median(gamma_excess)
            ),
            "fraction_positive": float(np.mean(gamma_excess > 0)),
        },
    }

    match_bad = (
        aggregate["median_passive_Rinput_match_factor"] > 1.50
        or aggregate["median_passive_soma_transfer_match_factor"] > 1.50
    )

    human_order_present = (
        aggregate["human"]["fraction_compact_order_above_5pct"] >= 0.667
    )
    compact_specific = (
        aggregate["human"][
            "median_compact_minus_dispersed_order_magnitude_log"
        ] >= threshold_log
        and aggregate["human"][
            "fraction_compact_more_order_sensitive_than_dispersed"
        ] >= 0.667
    )
    gamma_amplifies = (
        aggregate["gamma_specific"]["median_order_locality_excess_log"]
        >= threshold_log
        and aggregate["gamma_specific"]["fraction_positive"] >= 0.667
    )

    if match_bad:
        classification = "TIMING_ORDER_MATCH_INADEQUATE"
        interpretation = (
            "The timing-order contrast is not interpretable because the "
            "clustered/dispersed passive match is too loose."
        )
    elif not human_order_present:
        classification = "NO_ROBUST_COMPACT_TEMPORAL_ORDER_EFFECT"
        interpretation = (
            "Reversing the three-site arrival order does not change the HUMAN "
            "compact interaction by at least five percent on four of six "
            "branches after each order's own independent-superposition null."
        )
    elif not compact_specific:
        classification = "ORDER_EFFECT_NOT_COMPACT_SPECIFIC"
        interpretation = (
            "A temporal-order effect exists, but compact same-branch inputs are "
            "not robustly more order-sensitive than their mapped dispersed "
            "controls. The scaffold-specific order claim therefore fails."
        )
    elif gamma_amplifies:
        classification = "HUMAN_GAMMA_AMPLIFIES_COMPACT_TEMPORAL_ORDER"
        interpretation = (
            "The same compact dendritic sites respond differently to proximal- "
            "versus distal-first timing after linear superposition is removed, "
            "the effect is larger than in mapped dispersed controls, and the "
            "human NMDA gamma adds at least another five percent of compact-"
            "specific order sensitivity on the preregistered aggregate ruler."
        )
    else:
        classification = "COMPACT_TEMPORAL_ORDER_DEGREE_PRESENT"
        interpretation = (
            "The same compact dendritic sites retain a nonlinear temporal-order "
            "degree of freedom beyond their own linear time-shifted responses "
            "and beyond mapped dispersed controls. The extra order sensitivity "
            "is not robustly attributable to the human-vs-restmatched gamma "
            "difference alone."
        )

    summary = {
        "gate": 19,
        "object": "spatial scaffold x temporal order interaction",
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "same_pinned_released_cell1125": True,
            "cluster_sites_sorted_by_section_x": True,
            "proximal_to_distal_event_times_ms": forward_times.tolist(),
            "distal_to_proximal_event_times_ms": reverse_times.tolist(),
            "inter_site_lag_ms": float(args.lag_ms),
            "lag_locked_before_run": True,
            "primary_condition": "human",
            "attacker": "gamma062_restmatched",
            "independent_prediction": (
                "for every order and condition, sum three single-site traces "
                "run at those exact event times"
            ),
            "order_effect_threshold": (
                "absolute log forward/reverse interaction >= log(1.05)"
            ),
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "branches": result_rows,
        "stopping_line": (
            "A positive compact-specific order result earns a multi-pulse "
            "sequence discrimination test. It still does not earn growth. "
            "Failure closes the claim that dendritic path order and event "
            "order form a useful local computational degree of freedom here."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori Gate 19 — scaffold x temporal order")
    print()
    print(f"branches:                              {aggregate['branches']}")
    print(f"lag:                                   {aggregate['lag_ms']:.1f} ms")
    print(f"median passive Rinput match:           {aggregate['median_passive_Rinput_match_factor']:.3f}x")
    print(f"median passive soma-transfer match:    {aggregate['median_passive_soma_transfer_match_factor']:.3f}x")
    print()
    print("HUMAN")
    print(f"  compact order magnitude:             {aggregate['human']['median_compact_order_magnitude_fraction']:.4f}")
    print(f"  dispersed order magnitude:           {aggregate['human']['median_dispersed_order_magnitude_fraction']:.4f}")
    print(f"  compact branches >5%:                {aggregate['human']['fraction_compact_order_above_5pct']:.3f}")
    print(f"  compact-dispersed excess log:        {aggregate['human']['median_compact_minus_dispersed_order_magnitude_log']:.4f}")
    print(f"  compact more sensitive fraction:     {aggregate['human']['fraction_compact_more_order_sensitive_than_dispersed']:.3f}")
    print(f"  proximal->distal preferred:          {aggregate['human']['proximal_to_distal_preferred_fraction']:.3f}")
    print()
    print("rest-matched gamma=0.062")
    print(f"  compact order magnitude:             {aggregate['restmatched_gamma062']['median_compact_order_magnitude_fraction']:.4f}")
    print(f"  dispersed order magnitude:           {aggregate['restmatched_gamma062']['median_dispersed_order_magnitude_fraction']:.4f}")
    print(f"  compact-dispersed excess log:        {aggregate['restmatched_gamma062']['median_compact_minus_dispersed_order_magnitude_log']:.4f}")
    print()
    print("gamma-specific")
    print(f"  median order-locality excess log:    {aggregate['gamma_specific']['median_order_locality_excess_log']:.4f}")
    print(f"  positive fraction:                   {aggregate['gamma_specific']['fraction_positive']:.3f}")
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(result_rows) >= 4
    assert np.all(np.isfinite(human_compact_mag))
    assert np.all(np.isfinite(human_excess))
    assert np.all(np.isfinite(gamma_excess))
    assert classification in {
        "TIMING_ORDER_MATCH_INADEQUATE",
        "NO_ROBUST_COMPACT_TEMPORAL_ORDER_EFFECT",
        "ORDER_EFFECT_NOT_COMPACT_SPECIFIC",
        "COMPACT_TEMPORAL_ORDER_DEGREE_PRESENT",
        "HUMAN_GAMMA_AMPLIFIES_COMPACT_TEMPORAL_ORDER",
    }


if __name__ == "__main__":
    main()
