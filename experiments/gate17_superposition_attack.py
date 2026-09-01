from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from gate16_fci_dynamic_locality import (
    BASE_WEIGHT_US,
    FCI_COMMIT,
    HUMAN_GAMMA,
    MODEL_REL,
    V_INIT_MV,
    compact_midpoint_rows,
    dendritic_rows,
    git_head,
    import_builder,
    passive_site_features,
    section_groups,
)
from operaattori.compartment_match import greedy_dispersed_match


def block(v_mV: float, gamma: float = HUMAN_GAMMA) -> float:
    return float(1.0 / (1.0 + math.exp(-gamma * v_mV) / 3.57))


B_ZERO_GAMMA = 1.0 / (1.0 + 1.0 / 3.57)
HUMAN_RATIO = 0.00131 / 0.00088


def settle_baselines(syn_df, rows: np.ndarray) -> dict[int, float]:
    from neuron import h

    h.dt = 0.025
    h.finitialize(V_INIT_MV)
    h.continuerun(55.0)
    return {
        int(i): float(syn_df.iloc[int(i)]["segments"].v)
        for i in rows
    }


def configure_active(
    syn_df,
    active_rows: np.ndarray,
    multiplicity: float,
    condition: str,
    baseline_by_row: dict[int, float],
) -> None:
    for i in active_rows:
        row = syn_df.iloc[int(i)]
        syn = row["exc_synapses"]
        nc = row["exc_netcons"]

        if condition == "human":
            syn.gamma = HUMAN_GAMMA
            syn.NMDA_ratio = HUMAN_RATIO
        elif condition == "frozen":
            # Exact per-site frozen-at-rest effective NMDA conductance:
            # gamma=0 makes the released Mg gate constant at B_ZERO_GAMMA.
            # Rescaling the raw NMDA ratio by B(rest)/B_ZERO_GAMMA makes the
            # effective conductance equal to HUMAN at this site's actual
            # settled pre-event voltage.
            syn.gamma = 0.0
            syn.NMDA_ratio = (
                HUMAN_RATIO
                * block(float(baseline_by_row[int(i)]))
                / B_ZERO_GAMMA
            )
        else:
            raise ValueError(condition)

        nc.weight[0] = BASE_WEIGHT_US * float(multiplicity)


def run_trace(
    cell,
    syn_df,
    record_rows: np.ndarray,
    active_rows: np.ndarray,
    multiplicity: float,
    condition: str,
    baseline_by_row: dict[int, float],
    *,
    event_ms: float = 60.0,
    tstop_ms: float = 160.0,
) -> dict:
    from neuron import h

    record_rows = np.asarray(record_rows, dtype=int)
    active_rows = np.asarray(active_rows, dtype=int)
    configure_active(
        syn_df,
        active_rows,
        multiplicity,
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
    h.finitialize(V_INIT_MV)
    h.fcurrent()
    for i in active_rows:
        syn_df.iloc[int(i)]["exc_netcons"].event(float(event_ms))
    h.continuerun(float(tstop_ms))

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    local = np.stack(
        [np.asarray(v, dtype=float) for v in local_vecs],
        axis=0,
    )

    if not np.all(np.isfinite(soma)) or not np.all(np.isfinite(local)):
        raise FloatingPointError("non-finite Gate-17 trace")

    pre = (t >= event_ms - 10.0) & (t < event_ms - 1.0)
    post = (t >= event_ms) & (t <= event_ms + 90.0)
    soma_base = float(np.median(soma[pre]))
    local_base = np.median(local[:, pre], axis=1)

    return {
        "t": t[post],
        "soma_depol": soma[post] - soma_base,
        "local_depol": local[:, post] - local_base[:, None],
    }


def positive_auc(trace: np.ndarray, t: np.ndarray) -> float:
    return float(np.trapezoid(np.maximum(trace, 0.0), t))


def interaction_metrics(
    cell,
    syn_df,
    sites: np.ndarray,
    multiplicity: float,
    condition: str,
    baseline_by_row: dict[int, float],
) -> dict:
    sites = np.asarray(sites, dtype=int)

    together = run_trace(
        cell,
        syn_df,
        sites,
        sites,
        multiplicity,
        condition,
        baseline_by_row,
    )

    single_traces = [
        run_trace(
            cell,
            syn_df,
            sites,
            np.asarray([site], dtype=int),
            multiplicity,
            condition,
            baseline_by_row,
        )
        for site in sites
    ]

    t = together["t"]
    for one in single_traces:
        if not np.allclose(one["t"], t, rtol=0, atol=1e-12):
            raise RuntimeError("single-site trace grids differ")

    predicted_soma = np.sum(
        np.stack([x["soma_depol"] for x in single_traces], axis=0),
        axis=0,
    )
    predicted_local = np.sum(
        np.stack([x["local_depol"] for x in single_traces], axis=0),
        axis=0,
    )

    actual_local_mean = np.mean(together["local_depol"], axis=0)
    predicted_local_mean = np.mean(predicted_local, axis=0)

    actual_local_auc = positive_auc(actual_local_mean, t)
    predicted_local_auc = positive_auc(predicted_local_mean, t)
    actual_soma_auc = positive_auc(together["soma_depol"], t)
    predicted_soma_auc = positive_auc(predicted_soma, t)

    actual_peak = float(np.max(actual_local_mean))
    predicted_peak = float(np.max(predicted_local_mean))

    return {
        "actual_local_mean_auc_mV_ms": actual_local_auc,
        "independent_sum_local_mean_auc_mV_ms": predicted_local_auc,
        "local_interaction_ratio": float(
            actual_local_auc / (predicted_local_auc + 1e-30)
        ),
        "actual_soma_auc_mV_ms": actual_soma_auc,
        "independent_sum_soma_auc_mV_ms": predicted_soma_auc,
        "soma_interaction_ratio": float(
            actual_soma_auc / (predicted_soma_auc + 1e-30)
        ),
        "actual_local_mean_peak_mV": actual_peak,
        "independent_sum_local_mean_peak_mV": predicted_peak,
        "peak_interaction_ratio": float(
            actual_peak / (predicted_peak + 1e-30)
        ),
    }


def ratio(a: float, b: float) -> float:
    return float(a / (b + 1e-30))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=4)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument("--multiplicity", type=float, default=16.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/gate17/gate17_superposition_attack.json"),
    )
    args = ap.parse_args()

    if args.sites != 3:
        raise ValueError("Gate 17 protocol is locked to three sites")
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

    out_rows = []
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
            for condition in ("frozen", "human"):
                conditions[arrangement][condition] = interaction_metrics(
                    cell,
                    syn_df,
                    sites,
                    float(args.multiplicity),
                    condition,
                    baseline_by_row,
                )

        ch = conditions["clustered"]["human"]
        cf = conditions["clustered"]["frozen"]
        dh = conditions["dispersed"]["human"]
        df = conditions["dispersed"]["frozen"]

        cluster_nmda_interaction_gain = ratio(
            ch["local_interaction_ratio"],
            cf["local_interaction_ratio"],
        )
        dispersed_nmda_interaction_gain = ratio(
            dh["local_interaction_ratio"],
            df["local_interaction_ratio"],
        )
        locality = ratio(
            cluster_nmda_interaction_gain,
            dispersed_nmda_interaction_gain,
        )

        cluster_peak_gain = ratio(
            ch["peak_interaction_ratio"],
            cf["peak_interaction_ratio"],
        )
        dispersed_peak_gain = ratio(
            dh["peak_interaction_ratio"],
            df["peak_interaction_ratio"],
        )

        out_row = {
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
            "cluster_nmda_interaction_gain": cluster_nmda_interaction_gain,
            "dispersed_nmda_interaction_gain": dispersed_nmda_interaction_gain,
            "nmda_interaction_locality_index": locality,
            "nmda_peak_interaction_locality_index": ratio(
                cluster_peak_gain,
                dispersed_peak_gain,
            ),
            "conditions": conditions,
        }
        out_rows.append(out_row)

        print(
            f"[{bi+1:02d}/{len(candidates):02d}] {name} "
            f"span={span:.1f}um "
            f"matchR={match['median_z_ratio_factor']:.3f}x "
            f"matchT={match['median_transfer_ratio_factor']:.3f}x "
            f"cluster H/F-int={cluster_nmda_interaction_gain:.4f} "
            f"spread H/F-int={dispersed_nmda_interaction_gain:.4f} "
            f"L={locality:.4f}"
        )

    locality = np.asarray(
        [r["nmda_interaction_locality_index"] for r in out_rows],
        dtype=float,
    )
    peak_locality = np.asarray(
        [r["nmda_peak_interaction_locality_index"] for r in out_rows],
        dtype=float,
    )
    cluster_gain = np.asarray(
        [r["cluster_nmda_interaction_gain"] for r in out_rows],
        dtype=float,
    )
    dispersed_gain = np.asarray(
        [r["dispersed_nmda_interaction_gain"] for r in out_rows],
        dtype=float,
    )
    match_r = np.asarray(
        [r["passive_match"]["median_z_ratio_factor"] for r in out_rows],
        dtype=float,
    )
    match_t = np.asarray(
        [r["passive_match"]["median_transfer_ratio_factor"] for r in out_rows],
        dtype=float,
    )

    aggregate = {
        "branches": int(len(out_rows)),
        "virtual_synapses_per_site": float(args.multiplicity),
        "total_simultaneous_virtual_synapses": int(
            round(args.multiplicity * args.sites)
        ),
        "median_passive_Rinput_match_factor": float(np.median(match_r)),
        "median_passive_soma_transfer_match_factor": float(np.median(match_t)),
        "median_cluster_nmda_interaction_gain": float(np.median(cluster_gain)),
        "median_dispersed_nmda_interaction_gain": float(np.median(dispersed_gain)),
        "median_nmda_interaction_locality_index": float(np.median(locality)),
        "fraction_locality_over_1p05": float(np.mean(locality > 1.05)),
        "median_nmda_peak_interaction_locality_index": float(
            np.median(peak_locality)
        ),
        "settled_voltage_range_mV": [
            float(
                min(
                    min(r["settled_cluster_voltage_mV"])
                    for r in out_rows
                )
            ),
            float(
                max(
                    max(r["settled_cluster_voltage_mV"])
                    for r in out_rows
                )
            ),
        ],
    }

    if (
        aggregate["median_passive_Rinput_match_factor"] > 1.50
        or aggregate["median_passive_soma_transfer_match_factor"] > 1.50
    ):
        classification = "SUPERPOSITION_ATTACK_MATCH_INADEQUATE"
        interpretation = "Passive matched-site quality is insufficient."
    elif (
        aggregate["median_nmda_interaction_locality_index"] >= 1.05
        and aggregate["fraction_locality_over_1p05"] >= 0.50
    ):
        classification = "NMDA_INTERACTION_LOCALITY_SURVIVES_SUPERPOSITION_ATTACK"
        interpretation = (
            "The same-site simultaneous response departs from the sum of its own "
            "single-site traces more strongly under human NMDA on compact branches "
            "than in the matched dispersed control. Individual passive dynamics "
            "therefore do not explain the Gate-16 locality effect."
        )
    else:
        classification = "NO_NMDA_LOCALITY_AFTER_SUPERPOSITION_ATTACK"
        interpretation = (
            "The Gate-16 locality effect is not robust once each set is compared "
            "with the sum of its own independent single-site traces."
        )

    summary = {
        "gate": 17,
        "object": "same-site independent-superposition attack on Gate 16",
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "sites": int(args.sites),
            "cluster_span_um": float(args.cluster_span_um),
            "multiplicity_per_site": float(args.multiplicity),
            "frozen_block_uses_each_site_actual_settled_voltage": True,
            "independent_prediction": (
                "sum of three single-site voltage-depolarization time traces "
                "recorded at the same three readout sites"
            ),
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "branches": out_rows,
        "stopping_line": (
            "This attack removes individual passive temporal filtering as the "
            "explanation for Gate 16. Survival earns a timing-structure experiment; "
            "failure closes the locality claim."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print()
    print("Operaattori Gate 17 — same-site superposition attack")
    print()
    print(f"branches:                              {aggregate['branches']}")
    print(f"median passive Rinput match:           {aggregate['median_passive_Rinput_match_factor']:.3f}x")
    print(f"median passive soma-transfer match:    {aggregate['median_passive_soma_transfer_match_factor']:.3f}x")
    print(f"settled voltage range:                 {aggregate['settled_voltage_range_mV'][0]:.3f} .. {aggregate['settled_voltage_range_mV'][1]:.3f} mV")
    print(f"median cluster H/F interaction gain:   {aggregate['median_cluster_nmda_interaction_gain']:.4f}")
    print(f"median spread H/F interaction gain:    {aggregate['median_dispersed_nmda_interaction_gain']:.4f}")
    print(f"median interaction locality:           {aggregate['median_nmda_interaction_locality_index']:.4f}")
    print(f"fraction locality >1.05:               {aggregate['fraction_locality_over_1p05']:.3f}")
    print(f"median peak-interaction locality:      {aggregate['median_nmda_peak_interaction_locality_index']:.4f}")
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(out_rows) >= 4
    assert np.all(np.isfinite(locality))
    assert classification in {
        "SUPERPOSITION_ATTACK_MATCH_INADEQUATE",
        "NMDA_INTERACTION_LOCALITY_SURVIVES_SUPERPOSITION_ATTACK",
        "NO_NMDA_LOCALITY_AFTER_SUPERPOSITION_ATTACK",
    }


if __name__ == "__main__":
    main()
