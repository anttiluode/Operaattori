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
    compact_midpoint_rows,
    dendritic_rows,
    git_head,
    import_builder,
    interaction_metrics,
    passive_site_features,
    positive_auc,
    run_trace,
    section_groups,
    settle_baselines,
)


def cross_branch_metrics(
    cell,
    syn_df,
    sites_a: np.ndarray,
    sites_b: np.ndarray,
    multiplicity: float,
    condition: str,
    baseline_by_row: dict[int, float],
) -> dict:
    sites_a = np.asarray(sites_a, dtype=int)
    sites_b = np.asarray(sites_b, dtype=int)
    record_rows = np.concatenate([sites_a, sites_b])

    together = run_trace(
        cell,
        syn_df,
        record_rows,
        record_rows,
        multiplicity,
        condition,
        baseline_by_row,
    )
    a_only = run_trace(
        cell,
        syn_df,
        record_rows,
        sites_a,
        multiplicity,
        condition,
        baseline_by_row,
    )
    b_only = run_trace(
        cell,
        syn_df,
        record_rows,
        sites_b,
        multiplicity,
        condition,
        baseline_by_row,
    )

    t = together["t"]
    if (
        not np.allclose(a_only["t"], t, rtol=0, atol=1e-12)
        or not np.allclose(b_only["t"], t, rtol=0, atol=1e-12)
    ):
        raise RuntimeError("Gate-20 branch-pair time grids differ")

    predicted_local = a_only["local_depol"] + b_only["local_depol"]
    actual_local_mean = np.mean(together["local_depol"], axis=0)
    predicted_local_mean = np.mean(predicted_local, axis=0)

    actual_local_auc = positive_auc(actual_local_mean, t)
    predicted_local_auc = positive_auc(predicted_local_mean, t)
    interaction_ratio = float(
        actual_local_auc / (predicted_local_auc + 1e-30)
    )

    predicted_soma = a_only["soma_depol"] + b_only["soma_depol"]
    actual_soma_auc = positive_auc(together["soma_depol"], t)
    predicted_soma_auc = positive_auc(predicted_soma, t)

    above = together["soma_depol"] + float(
        # run_trace stores a baseline-subtracted soma trace, so recover the
        # approximate absolute voltage from a settled soma near -70 mV only
        # for a conservative spike guard below. We also use the much simpler
        # large-depolarization guard, which does not depend on exact baseline.
        -70.0
    ) >= -20.0
    spike_crossings = int(
        np.sum((~above[:-1]) & above[1:])
    ) if len(above) > 1 else 0
    large_soma_depol = bool(np.max(together["soma_depol"]) >= 40.0)

    return {
        "actual_local_mean_auc_mV_ms": float(actual_local_auc),
        "branch_alone_sum_local_mean_auc_mV_ms": float(
            predicted_local_auc
        ),
        "cross_branch_interaction_ratio": interaction_ratio,
        "cross_branch_interaction_magnitude_log": abs(
            float(math.log(interaction_ratio + 1e-30))
        ),
        "actual_soma_auc_mV_ms": float(actual_soma_auc),
        "branch_alone_sum_soma_auc_mV_ms": float(predicted_soma_auc),
        "soma_interaction_ratio": float(
            actual_soma_auc / (predicted_soma_auc + 1e-30)
        ),
        "conservative_soma_spike_crossings": spike_crossings,
        "large_soma_depolarization_guard": large_soma_depol,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=6)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument("--multiplicity", type=float, default=8.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/gate20/gate20_compartment_modularity.json"),
    )
    args = ap.parse_args()

    if args.sites != 3:
        raise ValueError("Gate 20 protocol is locked to three sites per branch")
    if abs(float(args.multiplicity) - 8.0) > 1e-12:
        raise ValueError("Gate 20 multiplicity is locked to 8 per site")
    if git_head(args.fci_root.resolve()) != FCI_COMMIT:
        raise RuntimeError("FCI checkout is not pinned")

    builder = import_builder(args.fci_root.resolve())
    cell, syn_df = builder.create_cell(
        path=str(args.fci_root.resolve() / MODEL_REL) + "/"
    )

    rows = dendritic_rows(syn_df)
    groups = section_groups(syn_df, rows)
    # The call settles the model and gives us the same deterministic section
    # ordering used in Gates 16-19. No nonlinear result is used for selection.
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
    candidates = candidates[: min(args.branches, len(candidates))]
    if len(candidates) != 6:
        raise RuntimeError("Gate 20 is locked to the six Gate-16 long sections")

    conditions = ("human", "gamma062_restmatched")
    within = {condition: [] for condition in conditions}

    for condition in conditions:
        for bi, branch in enumerate(candidates):
            metrics = interaction_metrics(
                cell,
                syn_df,
                branch["sites"],
                float(args.multiplicity),
                condition,
                baseline_by_row,
            )
            ratio = float(metrics["local_interaction_ratio"])
            row = {
                "branch_index": int(bi),
                "section": branch["section"],
                "span_um": branch["span_um"],
                "interaction_ratio": ratio,
                "interaction_magnitude_log": abs(
                    float(math.log(ratio + 1e-30))
                ),
                "metrics": metrics,
            }
            within[condition].append(row)
            print(
                f"within {condition:22s} [{bi+1}/6] "
                f"{branch['section']} I={ratio:.4f} "
                f"|logI|={row['interaction_magnitude_log']:.4f}"
            )

    pairs = list(itertools.combinations(range(len(candidates)), 2))
    pair_rows = []

    for pi, (ia, ib) in enumerate(pairs):
        branch_a = candidates[ia]
        branch_b = candidates[ib]
        row = {
            "pair_index": int(pi),
            "branch_a": int(ia),
            "branch_b": int(ib),
            "section_a": branch_a["section"],
            "section_b": branch_b["section"],
            "conditions": {},
        }

        for condition in conditions:
            cross = cross_branch_metrics(
                cell,
                syn_df,
                branch_a["sites"],
                branch_b["sites"],
                float(args.multiplicity),
                condition,
                baseline_by_row,
            )
            within_ref = 0.5 * (
                within[condition][ia]["interaction_magnitude_log"]
                + within[condition][ib]["interaction_magnitude_log"]
            )
            margin = float(
                within_ref
                - cross["cross_branch_interaction_magnitude_log"]
            )
            row["conditions"][condition] = {
                "within_reference_magnitude_log": float(within_ref),
                "cross": cross,
                "modularity_margin_log": margin,
            }

        row["human_minus_restmatched_modularity_margin_log"] = float(
            row["conditions"]["human"]["modularity_margin_log"]
            - row["conditions"]["gamma062_restmatched"][
                "modularity_margin_log"
            ]
        )
        pair_rows.append(row)

        h = row["conditions"]["human"]
        r = row["conditions"]["gamma062_restmatched"]
        print(
            f"pair [{pi+1:02d}/{len(pairs):02d}] "
            f"{ia}-{ib} "
            f"H within={h['within_reference_magnitude_log']:.4f} "
            f"cross={h['cross']['cross_branch_interaction_magnitude_log']:.4f} "
            f"margin={h['modularity_margin_log']:.4f} "
            f"R margin={r['modularity_margin_log']:.4f}"
        )

    threshold_log = float(math.log(1.05))

    human_within = np.asarray(
        [x["interaction_magnitude_log"] for x in within["human"]],
        dtype=float,
    )
    rest_within = np.asarray(
        [
            x["interaction_magnitude_log"]
            for x in within["gamma062_restmatched"]
        ],
        dtype=float,
    )
    human_cross = np.asarray(
        [
            p["conditions"]["human"]["cross"][
                "cross_branch_interaction_magnitude_log"
            ]
            for p in pair_rows
        ],
        dtype=float,
    )
    rest_cross = np.asarray(
        [
            p["conditions"]["gamma062_restmatched"]["cross"][
                "cross_branch_interaction_magnitude_log"
            ]
            for p in pair_rows
        ],
        dtype=float,
    )
    human_margin = np.asarray(
        [
            p["conditions"]["human"]["modularity_margin_log"]
            for p in pair_rows
        ],
        dtype=float,
    )
    rest_margin = np.asarray(
        [
            p["conditions"]["gamma062_restmatched"][
                "modularity_margin_log"
            ]
            for p in pair_rows
        ],
        dtype=float,
    )
    gamma_margin = human_margin - rest_margin
    spike_guard = np.asarray(
        [
            (
                p["conditions"]["human"]["cross"][
                    "conservative_soma_spike_crossings"
                ] > 0
                or p["conditions"]["human"]["cross"][
                    "large_soma_depolarization_guard"
                ]
            )
            for p in pair_rows
        ],
        dtype=bool,
    )

    aggregate = {
        "branches": 6,
        "branch_pairs": int(len(pair_rows)),
        "sites_per_branch": int(args.sites),
        "multiplicity_per_site": float(args.multiplicity),
        "synapses_per_active_branch": int(
            round(args.multiplicity * args.sites)
        ),
        "synapses_per_branch_pair": int(
            round(2.0 * args.multiplicity * args.sites)
        ),
        "human": {
            "median_within_branch_interaction_magnitude_log": float(
                np.median(human_within)
            ),
            "median_within_branch_interaction_fraction": float(
                np.exp(np.median(human_within)) - 1.0
            ),
            "fraction_branches_within_above_5pct": float(
                np.mean(human_within >= threshold_log)
            ),
            "median_cross_branch_interaction_magnitude_log": float(
                np.median(human_cross)
            ),
            "median_cross_branch_interaction_fraction": float(
                np.exp(np.median(human_cross)) - 1.0
            ),
            "median_modularity_margin_log": float(
                np.median(human_margin)
            ),
            "fraction_pairs_positive_modularity_margin": float(
                np.mean(human_margin > 0)
            ),
            "fraction_pairs_margin_above_log1p05": float(
                np.mean(human_margin >= threshold_log)
            ),
        },
        "restmatched_gamma062": {
            "median_within_branch_interaction_magnitude_log": float(
                np.median(rest_within)
            ),
            "median_cross_branch_interaction_magnitude_log": float(
                np.median(rest_cross)
            ),
            "median_modularity_margin_log": float(
                np.median(rest_margin)
            ),
        },
        "gamma_specific": {
            "median_human_minus_restmatched_modularity_margin_log": float(
                np.median(gamma_margin)
            ),
            "fraction_positive": float(np.mean(gamma_margin > 0)),
        },
        "fraction_human_pair_runs_spike_guard": float(
            np.mean(spike_guard)
        ),
    }

    if aggregate["fraction_human_pair_runs_spike_guard"] > 0.0:
        classification = "COMPARTMENT_MODULARITY_SPIKING_CONFOUNDED"
        interpretation = (
            "At least one branch-pair run crossed the conservative somatic "
            "spike guard, so the local modularity comparison is contaminated "
            "by global regenerative activity."
        )
    elif (
        aggregate["human"]["median_modularity_margin_log"]
        >= threshold_log
        and aggregate["human"][
            "fraction_pairs_positive_modularity_margin"
        ] >= 0.80
        and aggregate["human"][
            "fraction_branches_within_above_5pct"
        ] >= 0.667
    ):
        if (
            aggregate["gamma_specific"][
                "median_human_minus_restmatched_modularity_margin_log"
            ] >= threshold_log
            and aggregate["gamma_specific"]["fraction_positive"] >= 0.80
        ):
            classification = (
                "HUMAN_GAMMA_STRENGTHENS_SEMI_INDEPENDENT_COMPARTMENTS"
            )
            interpretation = (
                "Nonadditivity is substantially stronger inside compact "
                "branches than between branch compartments, and the human "
                "NMDA gamma adds a further robust modularity margin over the "
                "rest-matched gamma=0.062 attacker."
            )
        else:
            classification = "SEMI_INDEPENDENT_NONLINEAR_COMPARTMENTS_PRESENT"
            interpretation = (
                "Nonadditivity is substantially stronger within compact "
                "branches than between different branch compartments after "
                "each cross-branch pair is compared with the sum of its own "
                "two branch-alone responses. The human gamma is not required "
                "to explain the full modularity margin."
            )
    else:
        classification = "NO_ROBUST_COMPARTMENT_MODULARITY"
        interpretation = (
            "Within-branch nonlinear interaction is not robustly separated "
            "from cross-branch interaction by the preregistered modularity "
            "ruler."
        )

    summary = {
        "gate": 20,
        "object": "within-branch versus cross-branch nonlinear modularity",
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "same_six_long_compact_sections_as_gates16_19": True,
            "within_null": (
                "sum of the same branch's three single-site traces"
            ),
            "cross_null": (
                "sum of the two complete branch-cluster-alone traces; thus "
                "within-branch nonlinearity is already present in the null"
            ),
            "conditions": list(conditions),
            "multiplicity_per_site": float(args.multiplicity),
            "multiplicity_locked_before_run": True,
            "all_15_branch_pairs": True,
            "modularity_margin": (
                "mean within-branch |log interaction| of the two branches "
                "minus cross-branch |log interaction|"
            ),
            "threshold": "log(1.05), reusing the five-percent ruler",
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "within_branches": within,
        "pairs": pair_rows,
        "stopping_line": (
            "A positive modularity result earns a branch-subunit computation "
            "assay: same global input budget, different distributions across "
            "the measured compartments. Failure closes the simple "
            "semi-independent-subunit abstraction for this scaffold."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori Gate 20 — nonlinear compartment modularity")
    print()
    print(f"branches / pairs:                     6 / {len(pair_rows)}")
    print(f"synapses per branch / pair:           {aggregate['synapses_per_active_branch']} / {aggregate['synapses_per_branch_pair']}")
    print()
    print("HUMAN")
    print(f"  median within magnitude:            {aggregate['human']['median_within_branch_interaction_fraction']:.4f}")
    print(f"  branches within >5%:                {aggregate['human']['fraction_branches_within_above_5pct']:.3f}")
    print(f"  median cross magnitude:             {aggregate['human']['median_cross_branch_interaction_fraction']:.4f}")
    print(f"  median modularity margin log:       {aggregate['human']['median_modularity_margin_log']:.4f}")
    print(f"  pairs positive margin:              {aggregate['human']['fraction_pairs_positive_modularity_margin']:.3f}")
    print(f"  pairs margin >=log1.05:             {aggregate['human']['fraction_pairs_margin_above_log1p05']:.3f}")
    print()
    print("rest-matched gamma=0.062")
    print(f"  median modularity margin log:       {aggregate['restmatched_gamma062']['median_modularity_margin_log']:.4f}")
    print()
    print("gamma-specific")
    print(f"  median extra modularity log:        {aggregate['gamma_specific']['median_human_minus_restmatched_modularity_margin_log']:.4f}")
    print(f"  positive pairs:                     {aggregate['gamma_specific']['fraction_positive']:.3f}")
    print()
    print(f"pair runs with spike guard:           {aggregate['fraction_human_pair_runs_spike_guard']:.3f}")
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(pair_rows) == 15
    assert np.all(np.isfinite(human_within))
    assert np.all(np.isfinite(human_cross))
    assert np.all(np.isfinite(human_margin))
    assert classification in {
        "COMPARTMENT_MODULARITY_SPIKING_CONFOUNDED",
        "HUMAN_GAMMA_STRENGTHENS_SEMI_INDEPENDENT_COMPARTMENTS",
        "SEMI_INDEPENDENT_NONLINEAR_COMPARTMENTS_PRESENT",
        "NO_ROBUST_COMPARTMENT_MODULARITY",
    }


if __name__ == "__main__":
    main()
