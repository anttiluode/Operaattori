from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(EXPERIMENTS) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS))

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
from gate17_superposition_attack import HUMAN_RATIO, settle_baselines


DT_MS = 0.025
EVENT_MS = 60.0
TSTOP_MS = 160.0
POST_MS = 90.0
IMPULSE_NA = 0.001
LENGTH_SCALES = (0.80, 1.20)
BRANCHES = 6
SITES = 3
MULTIPLICITY = 8.0


def canonical_section_name(name: str) -> str:
    value = str(name)
    return value.split("].", 1)[1] if "]." in value else value


def make_cell(fci_root: Path):
    builder = import_builder(fci_root)
    return builder.create_cell(path=str(fci_root / MODEL_REL) + "/")


def recover_branches(cell, syn_df, cluster_span_um: float):
    rows = dendritic_rows(syn_df)
    groups = section_groups(syn_df, rows)
    passive_site_features(cell, syn_df, rows)

    candidates = []
    for name, row_ids in groups.items():
        clustered, span = compact_midpoint_rows(
            syn_df, row_ids, SITES, cluster_span_um
        )
        if len(clustered) != SITES:
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
                "canonical_section": canonical_section_name(name),
                "section_length_um": float(sec.L),
                "sites": clustered,
                "span_um": float(span),
                "site_x": [
                    float(syn_df.iloc[int(i)]["segments"].x)
                    for i in clustered
                ],
            }
        )

    candidates.sort(reverse=True, key=lambda x: x["section_length_um"])
    candidates = candidates[:BRANCHES]
    if len(candidates) != BRANCHES:
        raise RuntimeError("could not recover Gate-20 six-branch basis")
    return rows, candidates


def configure_human(syn_df, sites: np.ndarray) -> None:
    for i in np.asarray(sites, dtype=int):
        row = syn_df.iloc[int(i)]
        row["exc_synapses"].gamma = HUMAN_GAMMA
        row["exc_synapses"].NMDA_ratio = HUMAN_RATIO
        row["exc_netcons"].weight[0] = BASE_WEIGHT_US * MULTIPLICITY


def cluster_trace(cell, syn_df, sites: np.ndarray) -> dict:
    from neuron import h

    sites = np.asarray(sites, dtype=int)
    configure_human(syn_df, sites)

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
    ampa_vecs = [
        h.Vector().record(
            syn_df.iloc[int(i)]["exc_synapses"]._ref_i_AMPA
        )
        for i in sites
    ]
    nmda_vecs = [
        h.Vector().record(
            syn_df.iloc[int(i)]["exc_synapses"]._ref_i_NMDA
        )
        for i in sites
    ]

    h.dt = DT_MS
    h.finitialize(V_INIT_MV)
    h.fcurrent()
    for i in sites:
        syn_df.iloc[int(i)]["exc_netcons"].event(EVENT_MS)
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    ampa = np.stack(
        [np.asarray(v, dtype=float) for v in ampa_vecs], axis=0
    )
    nmda = np.stack(
        [np.asarray(v, dtype=float) for v in nmda_vecs], axis=0
    )

    pre = (t >= EVENT_MS - 10.0) & (t < EVENT_MS - 1.0)
    post = (t >= EVENT_MS) & (t <= EVENT_MS + POST_MS + 1e-9)
    if not np.any(pre) or not np.any(post):
        raise RuntimeError("missing factorization trace window")

    base = float(np.median(soma[pre]))
    soma_dep = soma[post] - base
    inward = -(ampa[:, post] + nmda[:, post])

    if not (
        np.all(np.isfinite(soma_dep))
        and np.all(np.isfinite(inward))
    ):
        raise FloatingPointError("non-finite factorization trace")

    return {
        "t": t[post] - EVENT_MS,
        "soma_depol": soma_dep,
        "site_inward_current_nA": inward,
        "soma_peak_absolute_mV": float(np.max(soma[post])),
        "spike_guard": bool(np.max(soma[post]) >= -20.0),
    }


def impulse_kernel(cell, sec, x: float) -> np.ndarray:
    from neuron import h

    stim = h.IClamp(float(x), sec=sec)
    stim.delay = EVENT_MS
    stim.dur = DT_MS
    stim.amp = IMPULSE_NA

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)

    h.dt = DT_MS
    h.finitialize(V_INIT_MV)
    h.fcurrent()
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    pre = (t >= EVENT_MS - 10.0) & (t < EVENT_MS - 1.0)
    post = (t >= EVENT_MS) & (t <= EVENT_MS + POST_MS + 1e-9)
    base = float(np.median(soma[pre]))
    return (soma[post] - base) / IMPULSE_NA


def transport_predict(currents: np.ndarray, kernels: np.ndarray) -> np.ndarray:
    currents = np.asarray(currents, dtype=float)
    kernels = np.asarray(kernels, dtype=float)
    if currents.shape != kernels.shape:
        raise ValueError(
            f"current/kernel shape mismatch {currents.shape} != {kernels.shape}"
        )
    n = currents.shape[1]
    out = np.zeros(n, dtype=float)
    for site in range(currents.shape[0]):
        out += np.convolve(
            currents[site], kernels[site], mode="full"
        )[:n]
    return out


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def positive_auc(v: np.ndarray, t: np.ndarray) -> float:
    return float(np.trapezoid(np.maximum(v, 0.0), t))


def trace_metrics(actual: np.ndarray, pred: np.ndarray, t: np.ndarray) -> dict:
    actual = np.asarray(actual, dtype=float)
    pred = np.asarray(pred, dtype=float)
    nrmse = rms(pred - actual) / (rms(actual) + 1e-30)

    apeak = float(np.max(actual))
    ppeak = float(np.max(pred))
    peak_rel = abs(ppeak - apeak) / (abs(apeak) + 1e-30)

    aauc = positive_auc(actual, t)
    pauc = positive_auc(pred, t)
    auc_rel = abs(pauc - aauc) / (abs(aauc) + 1e-30)

    ac = actual - float(np.mean(actual))
    pc = pred - float(np.mean(pred))
    corr = float(
        np.dot(ac, pc)
        / (np.linalg.norm(ac) * np.linalg.norm(pc) + 1e-30)
    )
    return {
        "nrmse": float(nrmse),
        "peak_relative_error": float(peak_rel),
        "positive_auc_relative_error": float(auc_rel),
        "correlation": corr,
        "actual_peak_mV": apeak,
        "predicted_peak_mV": ppeak,
        "actual_positive_auc_mV_ms": aauc,
        "predicted_positive_auc_mV_ms": pauc,
    }


def current_drift(original: np.ndarray, other: np.ndarray) -> dict:
    original = np.asarray(original, dtype=float)
    other = np.asarray(other, dtype=float)
    per_site = []
    for a, b in zip(original, other):
        per_site.append(
            {
                "nrmse": rms(b - a) / (rms(a) + 1e-30),
                "charge_relative_change": abs(
                    float(np.sum(b) - np.sum(a))
                ) / (abs(float(np.sum(a))) + 1e-30),
            }
        )
    return {
        "median_site_nrmse": float(
            np.median([x["nrmse"] for x in per_site])
        ),
        "median_site_charge_relative_change": float(
            np.median(
                [x["charge_relative_change"] for x in per_site]
            )
        ),
        "sites": per_site,
    }


def check_branch_identity(syn_df, sites: np.ndarray, expected: str):
    got = {
        canonical_section_name(
            syn_df.iloc[int(i)]["segments"].sec.name()
        )
        for i in sites
    }
    if got != {expected}:
        raise RuntimeError(
            f"branch identity mismatch: expected {expected}, got {got}"
        )


def compact_original_receipt(branch: dict, trace: dict, metrics: dict) -> dict:
    return {
        "branch_index": int(branch["branch_index"]),
        "section": branch["canonical_section"],
        "section_length_um": float(branch["section_length_um"]),
        "sites": [int(x) for x in branch["sites"]],
        "site_x": [float(x) for x in branch["site_x"]],
        "soma_peak_absolute_mV": float(trace["soma_peak_absolute_mV"]),
        "spike_guard": bool(trace["spike_guard"]),
        "reconstruction_metrics": metrics,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/operator_factorization/operator_factorization.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source is not pinned")

    cell0, syn0 = make_cell(fci_root)
    rows0, branches = recover_branches(
        cell0, syn0, args.cluster_span_um
    )
    settle_baselines(syn0, rows0)

    original = {}
    original_receipts = []
    original_reconstruction = []

    for bi, branch in enumerate(branches):
        branch = dict(branch)
        branch["branch_index"] = int(bi)
        sites = np.asarray(branch["sites"], dtype=int)
        tr = cluster_trace(cell0, syn0, sites)
        kernels = np.stack(
            [
                impulse_kernel(
                    cell0,
                    syn0.iloc[int(i)]["segments"].sec,
                    float(syn0.iloc[int(i)]["segments"].x),
                )
                for i in sites
            ],
            axis=0,
        )
        pred = transport_predict(
            tr["site_inward_current_nA"], kernels
        )
        metrics = trace_metrics(
            tr["soma_depol"], pred, tr["t"]
        )
        original[branch["canonical_section"]] = {
            "branch": branch,
            "trace": tr,
            "kernels": kernels,
            "reconstruction_metrics": metrics,
        }
        original_receipts.append(
            compact_original_receipt(branch, tr, metrics)
        )
        original_reconstruction.append(metrics["nrmse"])
        print(
            f"original [{bi+1}/6] {branch['canonical_section']} "
            f"transport NRMSE={metrics['nrmse']:.4f}"
        )

    holdouts = []
    spike_flags = []

    for scale in LENGTH_SCALES:
        for bi, branch in enumerate(branches):
            cell, syn = make_cell(fci_root)
            rows = dendritic_rows(syn)
            sites = np.asarray(branch["sites"], dtype=int)
            check_branch_identity(
                syn, sites, branch["canonical_section"]
            )

            sec = syn.iloc[int(sites[0])]["segments"].sec
            old_length = float(sec.L)
            if abs(
                old_length - branch["section_length_um"]
            ) > 1e-6:
                raise RuntimeError("fresh model length differs")
            sec.L = old_length * float(scale)
            new_length = float(sec.L)

            settle_baselines(syn, rows)
            actual = cluster_trace(cell, syn, sites)
            kernels = np.stack(
                [
                    impulse_kernel(
                        cell,
                        syn.iloc[int(i)]["segments"].sec,
                        float(syn.iloc[int(i)]["segments"].x),
                    )
                    for i in sites
                ],
                axis=0,
            )

            base = original[branch["canonical_section"]]
            t = actual["t"]
            if not np.allclose(
                base["trace"]["t"], t, rtol=0, atol=1e-12
            ):
                raise RuntimeError("holdout time grid changed")

            frozen_soma = base["trace"]["soma_depol"]
            factorized = transport_predict(
                base["trace"]["site_inward_current_nA"],
                kernels,
            )
            oracle = transport_predict(
                actual["site_inward_current_nA"],
                kernels,
            )

            frozen_metrics = trace_metrics(
                actual["soma_depol"], frozen_soma, t
            )
            factor_metrics = trace_metrics(
                actual["soma_depol"], factorized, t
            )
            oracle_metrics = trace_metrics(
                actual["soma_depol"], oracle, t
            )
            drift = current_drift(
                base["trace"]["site_inward_current_nA"],
                actual["site_inward_current_nA"],
            )
            spike_flags.append(actual["spike_guard"])

            holdouts.append(
                {
                    "scale": float(scale),
                    "branch_index": int(bi),
                    "section": branch["canonical_section"],
                    "section_length_original_um": old_length,
                    "section_length_holdout_um": new_length,
                    "frozen_soma_attacker": frozen_metrics,
                    "factorized": factor_metrics,
                    "transport_oracle": oracle_metrics,
                    "local_current_drift": drift,
                    "actual_spike_guard": actual["spike_guard"],
                }
            )

            print(
                f"holdout scale={scale:.2f} [{bi+1}/6] "
                f"{branch['canonical_section']} "
                f"frozen={frozen_metrics['nrmse']:.4f} "
                f"factor={factor_metrics['nrmse']:.4f} "
                f"oracle={oracle_metrics['nrmse']:.4f} "
                f"Jdrift={drift['median_site_nrmse']:.4f}"
            )

    original_err = np.asarray(
        original_reconstruction, dtype=float
    )
    frozen_err = np.asarray(
        [x["frozen_soma_attacker"]["nrmse"] for x in holdouts],
        dtype=float,
    )
    factor_err = np.asarray(
        [x["factorized"]["nrmse"] for x in holdouts],
        dtype=float,
    )
    oracle_err = np.asarray(
        [x["transport_oracle"]["nrmse"] for x in holdouts],
        dtype=float,
    )
    current_err = np.asarray(
        [x["local_current_drift"]["median_site_nrmse"] for x in holdouts],
        dtype=float,
    )

    aggregate = {
        "branches": BRANCHES,
        "holdout_scales": list(LENGTH_SCALES),
        "holdout_cases": int(len(holdouts)),
        "sites_per_branch": SITES,
        "multiplicity_per_site": MULTIPLICITY,
        "virtual_synapses_per_branch": int(SITES * MULTIPLICITY),
        "impulse_nA": IMPULSE_NA,
        "dt_ms": DT_MS,
        "median_original_reconstruction_nrmse": float(
            np.median(original_err)
        ),
        "median_frozen_soma_attacker_nrmse": float(
            np.median(frozen_err)
        ),
        "median_factorized_holdout_nrmse": float(
            np.median(factor_err)
        ),
        "median_transport_oracle_nrmse": float(
            np.median(oracle_err)
        ),
        "factorized_to_frozen_median_nrmse_ratio": float(
            np.median(factor_err) / (np.median(frozen_err) + 1e-30)
        ),
        "fraction_factorized_beats_frozen": float(
            np.mean(factor_err < frozen_err)
        ),
        "median_local_current_waveform_drift_nrmse": float(
            np.median(current_err)
        ),
        "fraction_actual_holdouts_spike_guard": float(
            np.mean(spike_flags)
        ),
    }

    transport_ok = (
        aggregate["median_original_reconstruction_nrmse"] <= 0.10
        and aggregate["median_transport_oracle_nrmse"] <= 0.10
    )
    factor_ok = (
        aggregate["median_factorized_holdout_nrmse"] <= 0.15
        and aggregate[
            "factorized_to_frozen_median_nrmse_ratio"
        ] <= 0.80
        and aggregate["fraction_factorized_beats_frozen"] >= 8.0 / 12.0
        and aggregate["fraction_actual_holdouts_spike_guard"] == 0.0
    )

    if transport_ok and factor_ok:
        classification = (
            "TRANSPORT_X_LOCAL_NONLINEAR_OPERATOR_FACTORIZATION"
        )
        interpretation = (
            "A nonlinear site-current operator measured only in the original "
            "geometry composes with geometry-specific linear transport kernels "
            "to predict held-out soma traces substantially better than carrying "
            "the original soma trace across geometry."
        )
    elif transport_ok:
        classification = (
            "TRANSPORT_FACTORIZATION_VALID_LOCAL_OPERATOR_NOT_PORTABLE"
        )
        interpretation = (
            "The site-current to soma transport operator reconstructs actual "
            "responses when supplied the correct held-out currents, but freezing "
            "the original local nonlinear current operator does not predict "
            "held-out geometry accurately enough."
        )
    else:
        classification = "LINEAR_TRANSPORT_COMPOSITION_INADEQUATE"
        interpretation = (
            "Even with oracle held-out synaptic current waveforms, the measured "
            "linear site-to-soma impulse kernels do not reconstruct the soma "
            "trace closely enough for this operator factorization."
        )

    summary = {
        "object": (
            "held-out composition of original local nonlinear synaptic-current "
            "operator with geometry-specific linear site-to-soma transport"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "same_gate20_six_branch_basis": True,
            "local_operator_training_geometry": 1.0,
            "heldout_length_scales": list(LENGTH_SCALES),
            "diameter_changed": False,
            "topology_changed": False,
            "normalized_synapse_x_changed": False,
            "human_nmda_kinetics_changed": False,
            "site_current_interface": "-(i_AMPA+i_NMDA)",
            "transport_kernel": (
                "somatic response to +0.001 nA one-dt IClamp at each site, "
                "divided by 0.001"
            ),
            "no_heldout_soma_fit": True,
            "thresholds_locked_before_run": {
                "original_reconstruction_median_nrmse_max": 0.10,
                "transport_oracle_median_nrmse_max": 0.10,
                "factorized_median_holdout_nrmse_max": 0.15,
                "factorized_to_frozen_median_nrmse_ratio_max": 0.80,
                "factorized_beats_frozen_fraction_min": 8.0 / 12.0,
                "spike_guard_fraction_max": 0.0,
            },
        },
        "aggregate": aggregate,
        "original_reconstruction": original_receipts,
        "holdouts": holdouts,
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not align traces, change impulse amplitude, tune geometry "
            "scales, add fitted correction factors or select favorable branches "
            "after seeing the result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print()
    print("Operaattori operator factorization audit")
    print()
    print(
        "median original reconstruction NRMSE: "
        f"{aggregate['median_original_reconstruction_nrmse']:.4f}"
    )
    print(
        "median frozen-soma attacker NRMSE:     "
        f"{aggregate['median_frozen_soma_attacker_nrmse']:.4f}"
    )
    print(
        "median factorized holdout NRMSE:       "
        f"{aggregate['median_factorized_holdout_nrmse']:.4f}"
    )
    print(
        "median transport-oracle NRMSE:         "
        f"{aggregate['median_transport_oracle_nrmse']:.4f}"
    )
    print(
        "factorized / frozen median error:      "
        f"{aggregate['factorized_to_frozen_median_nrmse_ratio']:.4f}"
    )
    print(
        "factorized beats frozen:               "
        f"{aggregate['fraction_factorized_beats_frozen']:.3f}"
    )
    print(
        "median local-current waveform drift:   "
        f"{aggregate['median_local_current_waveform_drift_nrmse']:.4f}"
    )
    print(
        "holdout spike guard:                   "
        f"{aggregate['fraction_actual_holdouts_spike_guard']:.3f}"
    )
    print(f"classification: {classification}")

    assert len(holdouts) == 12
    assert np.all(np.isfinite(original_err))
    assert np.all(np.isfinite(frozen_err))
    assert np.all(np.isfinite(factor_err))
    assert np.all(np.isfinite(oracle_err))


if __name__ == "__main__":
    main()
