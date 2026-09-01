from __future__ import annotations

import argparse
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.compartment_match import greedy_dispersed_match


FCI_COMMIT = "75ad8b4d81a7f51bf888b30650c543592340db06"
MODEL_REL = Path(
    "simulating_neurons/neuron_models/human/eyal/"
    "Human_L23_PC_0603_11_937_Eyal_passive_dends_simple_soma"
)
MORPHOLOGY = "2013_03_06_cell11_1125_H41_06.asc"
V_INIT_MV = -76.0
V_FREEZE_MV = -70.0
BASE_WEIGHT_US = 0.00088
HUMAN_GAMMA = 0.078


def git_head(root: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        text=True,
    ).strip()


def import_builder(fci_root: Path):
    sys.path.insert(0, str(fci_root))
    path = fci_root / MODEL_REL / "get_standard_model.py"
    spec = importlib.util.spec_from_file_location("fci_cell1125_builder_gate16", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dendritic_rows(syn_df) -> np.ndarray:
    out = []
    for i, row in syn_df.iterrows():
        name = row["segments"].sec.name()
        if ("dend" in name or "apic" in name) and "axon" not in name:
            out.append(int(i))
    return np.asarray(out, dtype=int)


def section_groups(syn_df, rows: np.ndarray):
    groups: dict[str, list[int]] = {}
    for i in rows:
        seg = syn_df.iloc[int(i)]["segments"]
        groups.setdefault(seg.sec.name(), []).append(int(i))
    for name in groups:
        groups[name].sort(key=lambda i: float(syn_df.iloc[i]["segments"].x))
    return groups


def compact_midpoint_rows(
    syn_df,
    row_ids: list[int],
    count: int,
    max_span_um: float,
) -> tuple[np.ndarray, float]:
    if len(row_ids) < count:
        return np.empty(0, dtype=int), 0.0

    rows = np.asarray(row_ids, dtype=int)
    xs = np.asarray(
        [float(syn_df.iloc[i]["segments"].x) for i in rows],
        dtype=float,
    )
    sec = syn_df.iloc[int(rows[0])]["segments"].sec
    L = float(sec.L)

    best = None
    for start in range(0, len(rows) - count + 1):
        chosen = rows[start : start + count]
        chosen_x = xs[start : start + count]
        span = float((chosen_x[-1] - chosen_x[0]) * L)
        if span > max_span_um + 1e-12:
            continue
        center_error = abs(float(np.mean(chosen_x)) - 0.5)
        key = (center_error, span)
        if best is None or key < best[0]:
            best = (key, chosen.copy(), span)

    if best is None:
        return np.empty(0, dtype=int), 0.0
    return best[1], float(best[2])


def passive_site_features(cell, syn_df, rows: np.ndarray):
    from neuron import h

    h.dt = 0.025
    h.finitialize(V_INIT_MV)
    h.continuerun(50.0)

    imp = h.Impedance()
    imp.loc(0.5, sec=cell.soma[0])
    imp.compute(0.0)

    n = len(syn_df)
    zinput = np.zeros(n, dtype=float)
    transfer = np.zeros(n, dtype=float)
    section_id = np.full(n, -1, dtype=np.int64)

    name_to_id: dict[str, int] = {}
    for i in rows:
        seg = syn_df.iloc[int(i)]["segments"]
        name = seg.sec.name()
        if name not in name_to_id:
            name_to_id[name] = len(name_to_id)
        section_id[int(i)] = name_to_id[name]
        zinput[int(i)] = float(imp.input(float(seg.x), sec=seg.sec))
        transfer[int(i)] = float(imp.transfer(float(seg.x), sec=seg.sec))

    return zinput, transfer, section_id, name_to_id


def frozen_nmda_ratio(
    original_ratio: float,
    *,
    gamma: float = HUMAN_GAMMA,
    v_freeze_mV: float = V_FREEZE_MV,
) -> float:
    # Released mechanism:
    # B(v) = 1 / (1 + exp(-gamma*v) * 1/3.57)
    # Setting gamma=0 makes B constant B0. Rescale NMDA_ratio so the effective
    # NMDA conductance equals the human mechanism at v_freeze for every v.
    b_rest = 1.0 / (1.0 + math.exp(-gamma * v_freeze_mV) / 3.57)
    b_zero_gamma = 1.0 / (1.0 + 1.0 / 3.57)
    return float(original_ratio * b_rest / b_zero_gamma)


def configure_sites(syn_df, site_rows: np.ndarray, multiplicity: float, condition: str):
    for i in site_rows:
        row = syn_df.iloc[int(i)]
        syn = row["exc_synapses"]
        nc = row["exc_netcons"]

        original_ratio = 0.00131 / 0.00088
        if condition == "human":
            syn.gamma = HUMAN_GAMMA
            syn.NMDA_ratio = original_ratio
        elif condition == "frozen":
            syn.gamma = 0.0
            syn.NMDA_ratio = frozen_nmda_ratio(original_ratio)
        else:
            raise ValueError(condition)

        # NET_RECEIVE is linear in event weight. Multiplying the event weight
        # is exactly equivalent to that many simultaneous identical synapses
        # arriving at this point process.
        nc.weight[0] = BASE_WEIGHT_US * float(multiplicity)


def run_arrangement(
    cell,
    syn_df,
    site_rows: np.ndarray,
    multiplicity: float,
    condition: str,
    *,
    event_ms: float = 60.0,
    tstop_ms: float = 160.0,
) -> dict:
    from neuron import h

    site_rows = np.asarray(site_rows, dtype=int)
    configure_sites(syn_df, site_rows, multiplicity, condition)

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)

    local_vecs = []
    nmda_vecs = []
    ampa_vecs = []
    for i in site_rows:
        row = syn_df.iloc[int(i)]
        seg = row["segments"]
        syn = row["exc_synapses"]
        local_vecs.append(h.Vector().record(seg._ref_v))
        nmda_vecs.append(h.Vector().record(syn._ref_i_NMDA))
        ampa_vecs.append(h.Vector().record(syn._ref_i_AMPA))

    h.dt = 0.025
    h.finitialize(V_INIT_MV)
    h.fcurrent()

    for i in site_rows:
        syn_df.iloc[int(i)]["exc_netcons"].event(float(event_ms))

    h.continuerun(float(tstop_ms))

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    local = np.stack([np.asarray(v, dtype=float) for v in local_vecs], axis=0)
    i_nmda = np.stack([np.asarray(v, dtype=float) for v in nmda_vecs], axis=0)
    i_ampa = np.stack([np.asarray(v, dtype=float) for v in ampa_vecs], axis=0)

    if (
        not np.all(np.isfinite(t))
        or not np.all(np.isfinite(soma))
        or not np.all(np.isfinite(local))
        or not np.all(np.isfinite(i_nmda))
        or not np.all(np.isfinite(i_ampa))
    ):
        raise FloatingPointError("non-finite released-model trace")

    pre = (t >= event_ms - 10.0) & (t < event_ms - 1.0)
    post = (t >= event_ms) & (t <= min(tstop_ms, event_ms + 90.0))
    if not np.any(pre) or not np.any(post):
        raise RuntimeError("missing baseline/post windows")

    soma_base = float(np.median(soma[pre]))
    local_base = np.median(local[:, pre], axis=1)

    soma_dep = soma - soma_base
    local_dep = local - local_base[:, None]
    local_mean_dep = np.mean(local_dep, axis=0)
    local_max_dep = np.max(local_dep, axis=0)

    tp = t[post]
    soma_auc = float(np.trapezoid(np.maximum(soma_dep[post], 0.0), tp))
    local_mean_auc = float(
        np.trapezoid(np.maximum(local_mean_dep[post], 0.0), tp)
    )
    local_max_auc = float(
        np.trapezoid(np.maximum(local_max_dep[post], 0.0), tp)
    )
    nmda_charge = float(
        np.trapezoid(np.sum(np.abs(i_nmda[:, post]), axis=0), tp)
    )
    ampa_charge = float(
        np.trapezoid(np.sum(np.abs(i_ampa[:, post]), axis=0), tp)
    )

    post_soma = soma[post]
    above = post_soma >= -20.0
    spike_count = int(np.sum((~above[:-1]) & above[1:])) if len(above) > 1 else 0

    return {
        "soma_baseline_mV": soma_base,
        "median_local_baseline_mV": float(np.median(local_base)),
        "soma_peak_depolarization_mV": float(np.max(soma_dep[post])),
        "local_mean_peak_depolarization_mV": float(np.max(local_mean_dep[post])),
        "local_max_peak_depolarization_mV": float(np.max(local_max_dep[post])),
        "soma_positive_auc_mV_ms": soma_auc,
        "local_mean_positive_auc_mV_ms": local_mean_auc,
        "local_max_positive_auc_mV_ms": local_max_auc,
        "NMDA_abs_current_auc_nA_ms": nmda_charge,
        "AMPA_abs_current_auc_nA_ms": ampa_charge,
        "soma_spike_count_minus20mV": spike_count,
    }


def ratio(a: float, b: float) -> float:
    return float(a / (b + 1e-30))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=6)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/gate16/gate16_dynamic_locality.json"),
    )
    args = ap.parse_args()

    if args.sites != 3:
        raise ValueError("Gate 16 protocol is locked to three matched sites")

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source is not at the pinned commit")

    builder = import_builder(fci_root)
    model_folder = fci_root / MODEL_REL
    cell, syn_df = builder.create_cell(path=str(model_folder) + "/")

    rows = dendritic_rows(syn_df)
    groups = section_groups(syn_df, rows)
    zinput, soma_transfer, run_ids, name_to_id = passive_site_features(
        cell, syn_df, rows
    )

    candidates = []
    for name, row_ids in groups.items():
        chosen, span = compact_midpoint_rows(
            syn_df, row_ids, args.sites, args.cluster_span_um
        )
        if len(chosen) != args.sites:
            continue
        sec = syn_df.iloc[int(chosen[0])]["segments"].sec
        candidates.append((float(sec.L), name, chosen, span))
    candidates.sort(reverse=True, key=lambda x: x[0])
    candidates = candidates[: min(args.branches, len(candidates))]
    if len(candidates) < 4:
        raise RuntimeError("too few released-model sections with compact 3-site windows")

    multiplicities = np.asarray([1.0, 4.0, 16.0])
    result_rows = []

    for bi, (section_length, name, clustered, actual_span) in enumerate(candidates):
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

        cluster_sections = {
            syn_df.iloc[int(i)]["segments"].sec.name() for i in clustered
        }
        dispersed_sections = {
            syn_df.iloc[int(i)]["segments"].sec.name() for i in dispersed
        }
        if len(cluster_sections) != 1 or len(dispersed_sections) != args.sites:
            raise RuntimeError("cluster/dispersal protocol violated")

        doses = []
        for mult in multiplicities:
            solved = {}
            for arrangement, sites in (
                ("clustered", clustered),
                ("dispersed", dispersed),
            ):
                solved[arrangement] = {}
                for condition in ("frozen", "human"):
                    solved[arrangement][condition] = run_arrangement(
                        cell,
                        syn_df,
                        sites,
                        float(mult),
                        condition,
                    )

            ch = solved["clustered"]["human"]
            cf = solved["clustered"]["frozen"]
            dh = solved["dispersed"]["human"]
            df = solved["dispersed"]["frozen"]

            local_hf_c = ratio(
                ch["local_mean_positive_auc_mV_ms"],
                cf["local_mean_positive_auc_mV_ms"],
            )
            local_hf_d = ratio(
                dh["local_mean_positive_auc_mV_ms"],
                df["local_mean_positive_auc_mV_ms"],
            )
            soma_hf_c = ratio(
                ch["soma_positive_auc_mV_ms"],
                cf["soma_positive_auc_mV_ms"],
            )
            soma_hf_d = ratio(
                dh["soma_positive_auc_mV_ms"],
                df["soma_positive_auc_mV_ms"],
            )

            doses.append(
                {
                    "virtual_synapses_per_site": float(mult),
                    "total_simultaneous_virtual_synapses": int(
                        round(float(mult) * args.sites)
                    ),
                    "local_auc_human_over_frozen_clustered": local_hf_c,
                    "local_auc_human_over_frozen_dispersed": local_hf_d,
                    "local_auc_locality_index": ratio(local_hf_c, local_hf_d),
                    "soma_auc_human_over_frozen_clustered": soma_hf_c,
                    "soma_auc_human_over_frozen_dispersed": soma_hf_d,
                    "soma_auc_locality_index": ratio(soma_hf_c, soma_hf_d),
                    "clustered": solved["clustered"],
                    "dispersed": solved["dispersed"],
                }
            )

        row = {
            "branch_index": int(bi),
            "section": name,
            "section_length_um": float(section_length),
            "actual_cluster_span_um": float(actual_span),
            "clustered_rows": clustered.tolist(),
            "dispersed_rows": dispersed.tolist(),
            "distinct_dispersed_sections": int(len(dispersed_sections)),
            "passive_match": match,
            "median_clustered_Rinput_Mohm": float(np.median(zinput[clustered])),
            "median_dispersed_Rinput_Mohm": float(np.median(zinput[dispersed])),
            "median_clustered_soma_transfer_Mohm": float(
                np.median(soma_transfer[clustered])
            ),
            "median_dispersed_soma_transfer_Mohm": float(
                np.median(soma_transfer[dispersed])
            ),
            "doses": doses,
        }
        result_rows.append(row)

        print(
            f"[{bi+1:02d}/{len(candidates):02d}] {name} "
            f"L={section_length:.1f}um span={actual_span:.1f}um "
            f"matchR={match['median_z_ratio_factor']:.3f}x "
            f"matchT={match['median_transfer_ratio_factor']:.3f}x "
            f"dynamic-L48={doses[-1]['local_auc_locality_index']:.4f}"
        )

    high_local = np.asarray(
        [r["doses"][-1]["local_auc_locality_index"] for r in result_rows],
        dtype=float,
    )
    high_soma = np.asarray(
        [r["doses"][-1]["soma_auc_locality_index"] for r in result_rows],
        dtype=float,
    )
    match_r = np.asarray(
        [r["passive_match"]["median_z_ratio_factor"] for r in result_rows],
        dtype=float,
    )
    match_t = np.asarray(
        [r["passive_match"]["median_transfer_ratio_factor"] for r in result_rows],
        dtype=float,
    )
    high_spikes = np.asarray(
        [
            r["doses"][-1]["clustered"]["human"]["soma_spike_count_minus20mV"]
            + r["doses"][-1]["dispersed"]["human"]["soma_spike_count_minus20mV"]
            for r in result_rows
        ],
        dtype=float,
    )

    dose_summary = []
    for di, mult in enumerate(multiplicities):
        vals = np.asarray(
            [r["doses"][di]["local_auc_locality_index"] for r in result_rows],
            dtype=float,
        )
        dose_summary.append(
            {
                "virtual_synapses_per_site": float(mult),
                "total_synapses": int(round(mult * args.sites)),
                "median_local_auc_locality_index": float(np.median(vals)),
                "fraction_over_1p05": float(np.mean(vals > 1.05)),
            }
        )

    aggregate = {
        "branches": int(len(result_rows)),
        "sites_per_arrangement": int(args.sites),
        "requested_cluster_span_um": float(args.cluster_span_um),
        "median_actual_cluster_span_um": float(
            np.median([r["actual_cluster_span_um"] for r in result_rows])
        ),
        "median_passive_Rinput_match_factor": float(np.median(match_r)),
        "median_passive_soma_transfer_match_factor": float(np.median(match_t)),
        "median_high_dose_local_auc_locality_index": float(np.median(high_local)),
        "fraction_high_dose_local_auc_locality_over_1p05": float(
            np.mean(high_local > 1.05)
        ),
        "median_high_dose_soma_auc_locality_index": float(np.median(high_soma)),
        "fraction_high_dose_conditions_with_somatic_spikes": float(
            np.mean(high_spikes > 0)
        ),
        "dose_summary": dose_summary,
    }

    if (
        aggregate["median_passive_Rinput_match_factor"] > 1.50
        or aggregate["median_passive_soma_transfer_match_factor"] > 1.50
    ):
        classification = "DYNAMIC_PASSIVE_MATCH_INADEQUATE"
        interpretation = (
            "The released-model temporal contrast is not interpretable because "
            "clustered and dispersed sites were not matched tightly enough."
        )
    elif aggregate["fraction_high_dose_conditions_with_somatic_spikes"] > 0.50:
        classification = "DYNAMIC_ASSAY_SPIKING_CONFOUNDED"
        interpretation = (
            "The high-dose temporal contrast is dominated by somatic spikes and "
            "cannot isolate the dendritic locality mechanism."
        )
    elif (
        aggregate["median_high_dose_local_auc_locality_index"] >= 1.05
        and aggregate["fraction_high_dose_local_auc_locality_over_1p05"] >= 0.60
    ):
        classification = "DYNAMIC_NMDA_LOCALITY_ADVANTAGE_PRESENT"
        interpretation = (
            "The released kinetic model restores an NMDA-specific same-branch "
            "locality advantage that was absent in the quasi-static reduction."
        )
    else:
        classification = "NO_DYNAMIC_NMDA_LOCALITY_ADVANTAGE"
        interpretation = (
            "Using the released cell-1125 kinetic synapses does not restore a "
            "robust NMDA-specific clustered-input advantage after passive matching."
        )

    summary = {
        "gate": 16,
        "object": "released time-domain matched NMDA locality audit",
        "fci_commit": FCI_COMMIT,
        "model": str(MODEL_REL),
        "morphology": MORPHOLOGY,
        "protocol": {
            "neuron_dt_ms": 0.025,
            "initial_voltage_mV": V_INIT_MV,
            "event_time_ms": 60.0,
            "analysis_window_ms": 90.0,
            "human_gamma_per_mV": HUMAN_GAMMA,
            "human_AMPA_event_weight_uS": BASE_WEIGHT_US,
            "frozen_control": (
                "released AMPANMDA_EMS with gamma=0 and NMDA_ratio rescaled "
                "so effective NMDA conductance equals human B(V=-70 mV)"
            ),
            "passive_matching": "NEURON Impedance input() + soma transfer() at 0 Hz",
            "not_FCI": True,
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "branches": result_rows,
        "stopping_line": (
            "Gate 16 is the last allowed rescue of the Gate-15 locality hypothesis. "
            "A negative released-model result closes this mechanistic line before "
            "growth; a positive result earns only a richer temporal perturbation."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print()
    print("Operaattori Gate 16 — released kinetic compartment locality")
    print()
    print(f"branches:                              {aggregate['branches']}")
    print(f"median compact span:                   {aggregate['median_actual_cluster_span_um']:.2f} um")
    print(f"median passive Rinput match:           {aggregate['median_passive_Rinput_match_factor']:.3f}x")
    print(f"median passive soma-transfer match:    {aggregate['median_passive_soma_transfer_match_factor']:.3f}x")
    print(f"median high-dose local-AUC locality:   {aggregate['median_high_dose_local_auc_locality_index']:.4f}")
    print(f"fraction local-AUC locality >1.05:     {aggregate['fraction_high_dose_local_auc_locality_over_1p05']:.3f}")
    print(f"median high-dose soma-AUC locality:    {aggregate['median_high_dose_soma_auc_locality_index']:.4f}")
    print(f"branches with high-dose soma spikes:   {aggregate['fraction_high_dose_conditions_with_somatic_spikes']:.3f}")
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(result_rows) >= 4
    assert np.all(np.isfinite(high_local))
    assert classification in {
        "DYNAMIC_PASSIVE_MATCH_INADEQUATE",
        "DYNAMIC_ASSAY_SPIKING_CONFOUNDED",
        "DYNAMIC_NMDA_LOCALITY_ADVANTAGE_PRESENT",
        "NO_DYNAMIC_NMDA_LOCALITY_ADVANTAGE",
    }


if __name__ == "__main__":
    main()
