from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"
AUDITS = ROOT / "audits"
for p in (ROOT, EXPERIMENTS, AUDITS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from gate16_fci_dynamic_locality import (
    FCI_COMMIT,
    dendritic_rows,
    git_head,
)
from gate17_superposition_attack import settle_baselines
from operator_factorization import (
    BRANCHES,
    IMPULSE_NA,
    LENGTH_SCALES,
    MULTIPLICITY,
    check_branch_identity,
    cluster_trace,
    current_drift,
    impulse_kernel,
    make_cell,
    noinput_soma_trace,
    recover_branches,
    trace_metrics,
    transport_predict,
)


PATTERNS = {
    "middle_single": (1,),
    "outer_pair": (0, 2),
    "triple": (0, 1, 2),
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/operator_factorization/cross_input_reuse.json"
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
    noinput0 = noinput_soma_trace(cell0)

    original = {}
    original_errors = []

    for bi, branch in enumerate(branches):
        all_sites = np.asarray(branch["sites"], dtype=int)
        kernels_all = np.stack(
            [
                impulse_kernel(
                    cell0,
                    syn0.iloc[int(i)]["segments"].sec,
                    float(syn0.iloc[int(i)]["segments"].x),
                    noinput0,
                )
                for i in all_sites
            ],
            axis=0,
        )

        original[branch["canonical_section"]] = {}
        for pattern_name, local_indices in PATTERNS.items():
            local_indices = np.asarray(local_indices, dtype=int)
            sites = all_sites[local_indices]
            tr = cluster_trace(cell0, syn0, sites, noinput0)
            kernels = kernels_all[local_indices]
            pred = transport_predict(
                tr["site_inward_current_nA"], kernels
            )
            metrics = trace_metrics(
                tr["soma_depol"], pred, tr["t"]
            )
            original[branch["canonical_section"]][pattern_name] = {
                "trace": tr,
                "metrics": metrics,
            }
            original_errors.append(metrics["nrmse"])
            print(
                f"original [{bi+1}/6] {branch['canonical_section']} "
                f"{pattern_name} NRMSE={metrics['nrmse']:.4f}"
            )

    holdouts = []
    spike_flags = []

    for scale in LENGTH_SCALES:
        for bi, branch in enumerate(branches):
            cell, syn = make_cell(fci_root)
            rows = dendritic_rows(syn)
            all_sites = np.asarray(branch["sites"], dtype=int)
            check_branch_identity(
                syn, all_sites, branch["canonical_section"]
            )

            sec = syn.iloc[int(all_sites[0])]["segments"].sec
            old_length = float(sec.L)
            if abs(old_length - branch["section_length_um"]) > 1e-6:
                raise RuntimeError("fresh model length differs")
            sec.L = old_length * float(scale)

            settle_baselines(syn, rows)
            noinput = noinput_soma_trace(cell)
            kernels_all = np.stack(
                [
                    impulse_kernel(
                        cell,
                        syn.iloc[int(i)]["segments"].sec,
                        float(syn.iloc[int(i)]["segments"].x),
                        noinput,
                    )
                    for i in all_sites
                ],
                axis=0,
            )

            for pattern_name, local_indices in PATTERNS.items():
                local_indices = np.asarray(local_indices, dtype=int)
                sites = all_sites[local_indices]
                kernels = kernels_all[local_indices]
                actual = cluster_trace(cell, syn, sites, noinput)
                base = original[
                    branch["canonical_section"]
                ][pattern_name]

                frozen = base["trace"]["soma_depol"]
                factorized = transport_predict(
                    base["trace"]["site_inward_current_nA"],
                    kernels,
                )
                oracle = transport_predict(
                    actual["site_inward_current_nA"],
                    kernels,
                )

                fm = trace_metrics(
                    actual["soma_depol"], frozen, actual["t"]
                )
                pm = trace_metrics(
                    actual["soma_depol"], factorized, actual["t"]
                )
                om = trace_metrics(
                    actual["soma_depol"], oracle, actual["t"]
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
                        "pattern": pattern_name,
                        "frozen_soma_attacker": fm,
                        "factorized": pm,
                        "transport_oracle": om,
                        "local_current_drift": drift,
                        "spike_guard": actual["spike_guard"],
                    }
                )
                print(
                    f"holdout {scale:.2f} [{bi+1}/6] "
                    f"{branch['canonical_section']} {pattern_name} "
                    f"frozen={fm['nrmse']:.4f} "
                    f"factor={pm['nrmse']:.4f} "
                    f"oracle={om['nrmse']:.4f}"
                )

    original_err = np.asarray(original_errors, dtype=float)
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

    pattern_medians = {}
    for pattern_name in PATTERNS:
        vals = np.asarray(
            [
                x["factorized"]["nrmse"]
                for x in holdouts
                if x["pattern"] == pattern_name
            ],
            dtype=float,
        )
        pattern_medians[pattern_name] = float(np.median(vals))

    aggregate = {
        "branches": BRANCHES,
        "patterns": list(PATTERNS),
        "holdout_scales": list(LENGTH_SCALES),
        "holdout_cases": int(len(holdouts)),
        "multiplicity_per_active_site": MULTIPLICITY,
        "impulse_nA": IMPULSE_NA,
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
        "pattern_factorized_median_nrmse": pattern_medians,
        "fraction_spike_guard": float(np.mean(spike_flags)),
    }

    pass_all = (
        aggregate["median_original_reconstruction_nrmse"] <= 0.10
        and aggregate["median_transport_oracle_nrmse"] <= 0.10
        and aggregate["median_factorized_holdout_nrmse"] <= 0.15
        and aggregate[
            "factorized_to_frozen_median_nrmse_ratio"
        ] <= 0.80
        and aggregate["fraction_factorized_beats_frozen"] >= 24.0 / 36.0
        and all(v <= 0.15 for v in pattern_medians.values())
        and aggregate["fraction_spike_guard"] == 0.0
    )

    if pass_all:
        classification = "TRANSPORT_OPERATOR_REUSES_ACROSS_INPUT_PATTERNS"
        interpretation = (
            "One geometry-specific site-to-soma transport operator composes "
            "with three distinct original-geometry local current operators per "
            "branch and predicts held-out metric perturbations better than "
            "carrying the original soma trace."
        )
    elif aggregate["median_transport_oracle_nrmse"] <= 0.10:
        classification = "TRANSPORT_REUSABLE_BUT_LOCAL_OPERATORS_GEOMETRY_SENSITIVE"
        interpretation = (
            "The transport module remains accurate across input patterns when "
            "supplied held-out currents, but freezing original local current "
            "operators is not robust enough across the full pattern panel."
        )
    else:
        classification = "TRANSPORT_OPERATOR_NOT_INPUT_INDEPENDENT_AT_REQUIRED_ACCURACY"
        interpretation = (
            "The shared transport kernels do not reconstruct held-out responses "
            "accurately across the locked input-pattern panel."
        )

    summary = {
        "object": (
            "reuse of one geometry-specific transport operator across three "
            "distinct local input-pattern operators"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "patterns": {
                k: list(v) for k, v in PATTERNS.items()
            },
            "holdout_scales": list(LENGTH_SCALES),
            "multiplicity_per_active_site": MULTIPLICITY,
            "shared_transport_per_geometry": True,
            "no_heldout_soma_fit": True,
            "thresholds_locked_before_run": {
                "original_median_nrmse_max": 0.10,
                "oracle_median_nrmse_max": 0.10,
                "factorized_median_nrmse_max": 0.15,
                "factorized_to_frozen_ratio_max": 0.80,
                "factorized_win_fraction_min": 24.0 / 36.0,
                "each_pattern_median_factorized_nrmse_max": 0.15,
                "spike_guard_fraction_max": 0.0,
            },
        },
        "aggregate": aggregate,
        "holdouts": holdouts,
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not introduce per-pattern transport gains, alignment, selected "
            "branches or alternative pattern panels after this result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori cross-input transport reuse audit")
    print()
    print(
        "median original reconstruction NRMSE: "
        f"{aggregate['median_original_reconstruction_nrmse']:.4f}"
    )
    print(
        "median frozen attacker NRMSE:          "
        f"{aggregate['median_frozen_soma_attacker_nrmse']:.4f}"
    )
    print(
        "median factorized NRMSE:               "
        f"{aggregate['median_factorized_holdout_nrmse']:.4f}"
    )
    print(
        "median transport oracle NRMSE:         "
        f"{aggregate['median_transport_oracle_nrmse']:.4f}"
    )
    print(
        "factorized/frozen:                     "
        f"{aggregate['factorized_to_frozen_median_nrmse_ratio']:.4f}"
    )
    print(
        "factorized beats frozen:               "
        f"{aggregate['fraction_factorized_beats_frozen']:.3f}"
    )
    for name, value in pattern_medians.items():
        print(f"pattern {name:16s} median:     {value:.4f}")
    print(
        "spike guard fraction:                  "
        f"{aggregate['fraction_spike_guard']:.3f}"
    )
    print(f"classification: {classification}")

    assert len(holdouts) == 36
    assert np.all(np.isfinite(factor_err))
    assert np.all(np.isfinite(oracle_err))


if __name__ == "__main__":
    main()
