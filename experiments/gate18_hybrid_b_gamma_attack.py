from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from gate17_superposition_attack import (
    FCI_COMMIT,
    HUMAN_GAMMA,
    HUMAN_RATIO,
    MODEL_REL,
    RAT_GAMMA,
    block,
    compact_midpoint_rows,
    dendritic_rows,
    git_head,
    import_builder,
    interaction_metrics,
    passive_site_features,
    section_groups,
    settle_baselines,
)
from operaattori.compartment_match import greedy_dispersed_match


def ratio(a: float, b: float) -> float:
    return float(a / (b + 1e-30))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=6)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument("--multiplicity", type=float, default=16.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/gate18/gate18_hybrid_b_gamma_attack.json"),
    )
    args = ap.parse_args()

    if args.sites != 3:
        raise ValueError("Gate 18 protocol is locked to three sites")
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
        sec = syn_df.iloc[int(clustered[0])]["segments"].sec
        candidates.append((float(sec.L), name, clustered, span))

    candidates.sort(reverse=True, key=lambda x: x[0])
    candidates = candidates[: min(args.branches, len(candidates))]
    if len(candidates) < 4:
        raise RuntimeError("too few compact sections")

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

        conditions = {}
        for arrangement, sites in (
            ("clustered", clustered),
            ("dispersed", dispersed),
        ):
            conditions[arrangement] = {}
            for condition in (
                "human",
                "hybrid_b",
                "gamma062_restmatched",
            ):
                conditions[arrangement][condition] = interaction_metrics(
                    cell,
                    syn_df,
                    sites,
                    float(args.multiplicity),
                    condition,
                    baseline_by_row,
                )

        ch = conditions["clustered"]["human"]["local_interaction_ratio"]
        cb = conditions["clustered"]["hybrid_b"]["local_interaction_ratio"]
        cr = conditions["clustered"]["gamma062_restmatched"][
            "local_interaction_ratio"
        ]

        dh = conditions["dispersed"]["human"]["local_interaction_ratio"]
        db = conditions["dispersed"]["hybrid_b"]["local_interaction_ratio"]
        dr = conditions["dispersed"]["gamma062_restmatched"][
            "local_interaction_ratio"
        ]

        exact_cluster_gain = ratio(ch, cb)
        exact_spread_gain = ratio(dh, db)
        exact_locality = ratio(exact_cluster_gain, exact_spread_gain)

        rest_cluster_gain = ratio(ch, cr)
        rest_spread_gain = ratio(dh, dr)
        rest_locality = ratio(rest_cluster_gain, rest_spread_gain)

        selected_rows = np.concatenate([clustered, dispersed])
        restmatch_scale = np.asarray(
            [
                block(float(baseline_by_row[int(i)]), HUMAN_GAMMA)
                / block(float(baseline_by_row[int(i)]), RAT_GAMMA)
                for i in selected_rows
            ],
            dtype=float,
        )

        result = {
            "branch_index": int(bi),
            "section": name,
            "section_length_um": float(section_length),
            "cluster_span_um": float(span),
            "clustered_rows": clustered.tolist(),
            "dispersed_rows": dispersed.tolist(),
            "passive_match": match,
            "settled_cluster_voltage_mV": [
                float(baseline_by_row[int(i)]) for i in clustered
            ],
            "settled_dispersed_voltage_mV": [
                float(baseline_by_row[int(i)]) for i in dispersed
            ],
            "restmatched_gamma062_nmda_ratio_scale_vs_human": {
                "median": float(np.median(restmatch_scale)),
                "min": float(np.min(restmatch_scale)),
                "max": float(np.max(restmatch_scale)),
            },
            "paper_hybrid_b": {
                "cluster_human_over_hybridB_interaction": exact_cluster_gain,
                "spread_human_over_hybridB_interaction": exact_spread_gain,
                "locality_index": exact_locality,
            },
            "restmatched_gamma062": {
                "cluster_human_over_gamma062_interaction": rest_cluster_gain,
                "spread_human_over_gamma062_interaction": rest_spread_gain,
                "locality_index": rest_locality,
            },
            "conditions": conditions,
        }
        result_rows.append(result)

        print(
            f"[{bi+1:02d}/{len(candidates):02d}] {name} "
            f"span={span:.1f}um "
            f"matchR={match['median_z_ratio_factor']:.3f}x "
            f"matchT={match['median_transfer_ratio_factor']:.3f}x "
            f"L-paperB={exact_locality:.4f} "
            f"L-restmatched={rest_locality:.4f}"
        )

    exact_locality = np.asarray(
        [r["paper_hybrid_b"]["locality_index"] for r in result_rows],
        dtype=float,
    )
    rest_locality = np.asarray(
        [r["restmatched_gamma062"]["locality_index"] for r in result_rows],
        dtype=float,
    )
    exact_cluster = np.asarray(
        [
            r["paper_hybrid_b"]["cluster_human_over_hybridB_interaction"]
            for r in result_rows
        ],
        dtype=float,
    )
    exact_spread = np.asarray(
        [
            r["paper_hybrid_b"]["spread_human_over_hybridB_interaction"]
            for r in result_rows
        ],
        dtype=float,
    )
    rest_cluster = np.asarray(
        [
            r["restmatched_gamma062"][
                "cluster_human_over_gamma062_interaction"
            ]
            for r in result_rows
        ],
        dtype=float,
    )
    rest_spread = np.asarray(
        [
            r["restmatched_gamma062"][
                "spread_human_over_gamma062_interaction"
            ]
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
    rest_scales = np.asarray(
        [
            r["restmatched_gamma062_nmda_ratio_scale_vs_human"]["median"]
            for r in result_rows
        ],
        dtype=float,
    )

    aggregate = {
        "branches": int(len(result_rows)),
        "virtual_synapses_per_site": float(args.multiplicity),
        "total_simultaneous_virtual_synapses": int(
            round(args.multiplicity * args.sites)
        ),
        "human_gamma_per_mV": HUMAN_GAMMA,
        "hybrid_b_gamma_per_mV": RAT_GAMMA,
        "human_raw_nmda_ratio": HUMAN_RATIO,
        "median_passive_Rinput_match_factor": float(np.median(match_r)),
        "median_passive_soma_transfer_match_factor": float(np.median(match_t)),
        "paper_hybrid_b": {
            "median_cluster_human_over_hybridB_interaction": float(
                np.median(exact_cluster)
            ),
            "median_spread_human_over_hybridB_interaction": float(
                np.median(exact_spread)
            ),
            "median_locality_index": float(np.median(exact_locality)),
            "fraction_locality_over_1p05": float(
                np.mean(exact_locality > 1.05)
            ),
        },
        "restmatched_gamma062": {
            "median_nmda_ratio_scale_vs_human": float(
                np.median(rest_scales)
            ),
            "median_cluster_human_over_gamma062_interaction": float(
                np.median(rest_cluster)
            ),
            "median_spread_human_over_gamma062_interaction": float(
                np.median(rest_spread)
            ),
            "median_locality_index": float(np.median(rest_locality)),
            "fraction_locality_over_1p05": float(
                np.mean(rest_locality > 1.05)
            ),
        },
    }

    match_bad = (
        aggregate["median_passive_Rinput_match_factor"] > 1.50
        or aggregate["median_passive_soma_transfer_match_factor"] > 1.50
    )

    exact_pass = (
        aggregate["paper_hybrid_b"]["median_locality_index"] >= 1.05
        and aggregate["paper_hybrid_b"]["fraction_locality_over_1p05"] >= 0.50
    )
    rest_pass = (
        aggregate["restmatched_gamma062"]["median_locality_index"] >= 1.05
        and aggregate["restmatched_gamma062"][
            "fraction_locality_over_1p05"
        ] >= 0.50
    )

    if match_bad:
        classification = "HYBRID_B_ATTACK_MATCH_INADEQUATE"
        interpretation = (
            "The Hybrid-B gamma attack cannot be interpreted because the "
            "clustered/dispersed passive match is too loose."
        )
    elif rest_pass:
        classification = (
            "HUMAN_GAMMA_LOCALITY_SURVIVES_RESTMATCHED_HYBRID_B"
        )
        interpretation = (
            "The compact-branch interaction advantage survives not only the "
            "paper's Hybrid B control but also a gamma=0.062 control whose "
            "effective NMDA strength is matched to HUMAN at each site's own "
            "resting voltage. The remaining difference is therefore tied to "
            "the voltage-dependence curve after depolarization, not merely "
            "to resting NMDA strength."
        )
    elif exact_pass:
        classification = (
            "PAPER_HYBRID_B_DIFFERS_BUT_RESTMATCHED_SLOPE_ATTACK_FAILS"
        )
        interpretation = (
            "HUMAN differs from the paper's raw Hybrid B in a locality-shaped "
            "way, but the effect disappears when gamma=0.062 is rest-matched "
            "site by site. Resting effective NMDA strength is therefore a "
            "sufficient explanation for this assay."
        )
    else:
        classification = "NO_ROBUST_HUMAN_GAMMA_LOCALITY_ADVANTAGE"
        interpretation = (
            "The Gate-17 interaction locality does not robustly distinguish "
            "human gamma=0.078 from the smaller gamma=0.062 controls."
        )

    summary = {
        "gate": 18,
        "object": (
            "paper Hybrid-B and rest-matched gamma-slope attack on Gate 17"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "sites": int(args.sites),
            "cluster_span_um": float(args.cluster_span_um),
            "multiplicity_per_site": float(args.multiplicity),
            "human_gamma": HUMAN_GAMMA,
            "paper_hybrid_b_gamma": RAT_GAMMA,
            "paper_hybrid_b_raw_nmda_ratio_unchanged": True,
            "restmatched_control": (
                "gamma=0.062 with each site's NMDA_ratio rescaled so its "
                "effective NMDA conductance at that site's actual settled "
                "pre-event voltage equals HUMAN"
            ),
            "interaction_ruler": (
                "simultaneous local-voltage AUC divided by sum of the same "
                "three sites' independent single-site traces"
            ),
            "threshold_reused_from_gate17": 1.05,
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "branches": result_rows,
        "stopping_line": (
            "A positive rest-matched result earns a timing/order perturbation "
            "within the same released model. It still does not earn growth. "
            "Failure closes the claim that the human gamma slope itself "
            "creates the compact interaction advantage."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori Gate 18 — Hybrid B / gamma-slope attack")
    print()
    print(f"branches:                              {aggregate['branches']}")
    print(f"median passive Rinput match:           {aggregate['median_passive_Rinput_match_factor']:.3f}x")
    print(f"median passive soma-transfer match:    {aggregate['median_passive_soma_transfer_match_factor']:.3f}x")
    print()
    print("paper Hybrid B (same raw human conductances, gamma=0.062)")
    print(f"  median cluster H/B interaction:      {aggregate['paper_hybrid_b']['median_cluster_human_over_hybridB_interaction']:.4f}")
    print(f"  median spread H/B interaction:       {aggregate['paper_hybrid_b']['median_spread_human_over_hybridB_interaction']:.4f}")
    print(f"  median locality:                     {aggregate['paper_hybrid_b']['median_locality_index']:.4f}")
    print(f"  fraction locality >1.05:             {aggregate['paper_hybrid_b']['fraction_locality_over_1p05']:.3f}")
    print()
    print("rest-matched gamma=0.062")
    print(f"  median NMDA ratio scale vs human:    {aggregate['restmatched_gamma062']['median_nmda_ratio_scale_vs_human']:.4f}")
    print(f"  median cluster H/R interaction:      {aggregate['restmatched_gamma062']['median_cluster_human_over_gamma062_interaction']:.4f}")
    print(f"  median spread H/R interaction:       {aggregate['restmatched_gamma062']['median_spread_human_over_gamma062_interaction']:.4f}")
    print(f"  median locality:                     {aggregate['restmatched_gamma062']['median_locality_index']:.4f}")
    print(f"  fraction locality >1.05:             {aggregate['restmatched_gamma062']['fraction_locality_over_1p05']:.3f}")
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(result_rows) >= 4
    assert np.all(np.isfinite(exact_locality))
    assert np.all(np.isfinite(rest_locality))
    assert classification in {
        "HYBRID_B_ATTACK_MATCH_INADEQUATE",
        "HUMAN_GAMMA_LOCALITY_SURVIVES_RESTMATCHED_HYBRID_B",
        "PAPER_HYBRID_B_DIFFERS_BUT_RESTMATCHED_SLOPE_ATTACK_FAILS",
        "NO_ROBUST_HUMAN_GAMMA_LOCALITY_ADVANTAGE",
    }


if __name__ == "__main__":
    main()
