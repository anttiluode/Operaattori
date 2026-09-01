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
    FCI_COMMIT,
    MODEL_REL,
    compact_midpoint_rows,
    dendritic_rows,
    git_head,
    import_builder,
    passive_site_features,
    section_groups,
)
from gate17_superposition_attack import (
    positive_auc,
    run_trace,
    settle_baselines,
)


LENGTH_SCALE = 1.20
MULTIPLICITY = 8.0
SITES = 3
BRANCHES = 6
THRESHOLD_LOG = math.log(1.05)


def recover_branches(cell, syn_df, rows, cluster_span_um: float):
    groups = section_groups(syn_df, rows)
    # Same deterministic passive settle/order step used by Gate 20.
    zinput, transfer, _run_ids, _name_to_id = passive_site_features(
        cell, syn_df, rows
    )

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
                "section_length_um": float(sec.L),
                "sites": clustered,
                "span_um": float(span),
                "median_zinput_Mohm": float(np.median(zinput[clustered])),
                "median_transfer_Mohm": float(np.median(transfer[clustered])),
            }
        )

    candidates.sort(reverse=True, key=lambda x: x["section_length_um"])
    candidates = candidates[:BRANCHES]
    if len(candidates) != BRANCHES:
        raise RuntimeError("could not recover Gate-20 six-branch basis")
    return candidates, zinput, transfer


def interaction_with_guard(
    cell,
    syn_df,
    sites: np.ndarray,
    baseline_by_row: dict[int, float],
):
    sites = np.asarray(sites, dtype=int)
    together = run_trace(
        cell,
        syn_df,
        sites,
        sites,
        MULTIPLICITY,
        "human",
        baseline_by_row,
    )
    singles = [
        run_trace(
            cell,
            syn_df,
            sites,
            np.asarray([site], dtype=int),
            MULTIPLICITY,
            "human",
            baseline_by_row,
        )
        for site in sites
    ]

    t = together["t"]
    for one in singles:
        if not np.allclose(one["t"], t, rtol=0, atol=1e-12):
            raise RuntimeError("metric audit trace grids differ")

    predicted_local = np.sum(
        np.stack([x["local_depol"] for x in singles], axis=0),
        axis=0,
    )
    actual_mean = np.mean(together["local_depol"], axis=0)
    predicted_mean = np.mean(predicted_local, axis=0)
    actual_auc = positive_auc(actual_mean, t)
    predicted_auc = positive_auc(predicted_mean, t)
    interaction = float(actual_auc / (predicted_auc + 1e-30))

    absolute_soma = (
        together["soma_depol"] + float(together["soma_baseline_mV"])
    )
    spike_guard = bool(np.max(absolute_soma) >= -20.0)

    return {
        "local_interaction_ratio": interaction,
        "actual_local_auc_mV_ms": float(actual_auc),
        "independent_local_auc_mV_ms": float(predicted_auc),
        "soma_peak_absolute_mV": float(np.max(absolute_soma)),
        "spike_guard": spike_guard,
    }


def relative_log_factor(a: float, b: float) -> float:
    return float(abs(math.log((b + 1e-30) / (a + 1e-30))))


def make_cell(fci_root: Path):
    builder = import_builder(fci_root)
    return builder.create_cell(path=str(fci_root / MODEL_REL) + "/")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/metric_nonlinearity/metric_nonlinearity.json"),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source is not pinned")

    # Untouched reference model: recover the exact branch identities and
    # measure all original nonlinear interactions once.
    cell0, syn0 = make_cell(fci_root)
    rows0 = dendritic_rows(syn0)
    branches, z0_all, transfer0_all = recover_branches(
        cell0, syn0, rows0, args.cluster_span_um
    )
    baseline0 = settle_baselines(syn0, rows0)

    original = {}
    for branch in branches:
        original[branch["section"]] = interaction_with_guard(
            cell0, syn0, branch["sites"], baseline0
        )

    rows_out = []
    for bi, branch in enumerate(branches):
        # Fresh model for each metric intervention prevents one branch edit
        # from contaminating the next branch.
        cell, syn = make_cell(fci_root)
        rows = dendritic_rows(syn)

        sites = np.asarray(branch["sites"], dtype=int)
        names = {
            syn.iloc[int(i)]["segments"].sec.name()
            for i in sites
        }
        if names != {branch["section"]}:
            raise RuntimeError(
                f"site identity changed for {branch['section']}: {names}"
            )

        sec = syn.iloc[int(sites[0])]["segments"].sec
        old_length = float(sec.L)
        if abs(old_length - branch["section_length_um"]) > 1e-6:
            raise RuntimeError("fresh model section length differs from reference")

        # Intrinsic metric edit. With pt3d geometry NEURON may rescale the
        # section's stored 3-D shape as a consequence, but the causal quantity
        # being manipulated here is sec.L. Diameter/topology and normalized
        # synapse x locations are not changed.
        sec.L = old_length * LENGTH_SCALE
        new_length = float(sec.L)
        if abs(new_length / old_length - LENGTH_SCALE) > 1e-9:
            raise RuntimeError("NEURON did not apply the locked 20% L change")

        z1_all, transfer1_all, _run_ids, _name_to_id = passive_site_features(
            cell, syn, rows
        )
        baseline1 = settle_baselines(syn, rows)
        stretched = interaction_with_guard(
            cell, syn, sites, baseline1
        )

        original_metrics = original[branch["section"]]
        effect_log = relative_log_factor(
            original_metrics["local_interaction_ratio"],
            stretched["local_interaction_ratio"],
        )

        z0 = float(np.median(z0_all[sites]))
        z1 = float(np.median(z1_all[sites]))
        t0 = float(np.median(transfer0_all[sites]))
        t1 = float(np.median(transfer1_all[sites]))
        passive_effect = max(
            relative_log_factor(z0, z1),
            relative_log_factor(t0, t1),
        )

        row = {
            "branch_index": int(bi),
            "section": branch["section"],
            "sites": sites.tolist(),
            "compact_span_original_um": branch["span_um"],
            "section_length_original_um": old_length,
            "section_length_stretched_um": new_length,
            "length_factor": float(new_length / old_length),
            "median_input_impedance_original_Mohm": z0,
            "median_input_impedance_stretched_Mohm": z1,
            "median_soma_transfer_original_Mohm": t0,
            "median_soma_transfer_stretched_Mohm": t1,
            "passive_site_effect_magnitude_log": passive_effect,
            "original": original_metrics,
            "stretched": stretched,
            "nonlinear_interaction_effect_magnitude_log": effect_log,
            "nonlinear_interaction_effect_factor": float(math.exp(effect_log)),
        }
        rows_out.append(row)

        print(
            f"[{bi+1}/6] {branch['section']} "
            f"L {old_length:.1f}->{new_length:.1f}um "
            f"passive={math.exp(passive_effect):.4f}x "
            f"I {original_metrics['local_interaction_ratio']:.4f}"
            f"->{stretched['local_interaction_ratio']:.4f} "
            f"effect={math.exp(effect_log):.4f}x"
        )

    nonlinear = np.asarray(
        [r["nonlinear_interaction_effect_magnitude_log"] for r in rows_out],
        dtype=float,
    )
    passive = np.asarray(
        [r["passive_site_effect_magnitude_log"] for r in rows_out],
        dtype=float,
    )
    spike_guard = np.asarray(
        [r["original"]["spike_guard"] or r["stretched"]["spike_guard"] for r in rows_out],
        dtype=bool,
    )

    aggregate = {
        "branches": BRANCHES,
        "sites_per_branch": SITES,
        "multiplicity_per_site": MULTIPLICITY,
        "virtual_synapses_per_branch": int(SITES * MULTIPLICITY),
        "length_scale": LENGTH_SCALE,
        "median_passive_site_effect_magnitude_log": float(np.median(passive)),
        "median_passive_site_effect_factor": float(math.exp(np.median(passive))),
        "fraction_passive_effect_over_1pct": float(
            np.mean(passive >= math.log(1.01))
        ),
        "median_nonlinear_interaction_effect_magnitude_log": float(
            np.median(nonlinear)
        ),
        "median_nonlinear_interaction_effect_factor": float(
            math.exp(np.median(nonlinear))
        ),
        "fraction_branches_nonlinear_effect_over_5pct": float(
            np.mean(nonlinear >= THRESHOLD_LOG)
        ),
        "fraction_runs_spike_guard": float(np.mean(spike_guard)),
    }

    passive_ok = (
        aggregate["median_passive_site_effect_magnitude_log"]
        >= math.log(1.01)
    )
    nonlinear_ok = (
        aggregate["median_nonlinear_interaction_effect_magnitude_log"]
        >= THRESHOLD_LOG
        and aggregate["fraction_branches_nonlinear_effect_over_5pct"]
        >= 4.0 / 6.0
        and aggregate["fraction_runs_spike_guard"] == 0.0
    )

    if passive_ok and nonlinear_ok:
        classification = "INTRINSIC_METRIC_MODULATES_NONLINEAR_BRANCH_INTERACTION"
        interpretation = (
            "The same 20% local cable-length edit that changes passive site "
            "properties also changes the exact nonlinear three-site interaction "
            "ratio across most of the established compact branches."
        )
    elif passive_ok:
        classification = "METRIC_CHANGES_PASSIVE_TRANSPORT_BUT_NOT_NONLINEAR_RATIO"
        interpretation = (
            "The local cable-length intervention is electrically real by the "
            "passive positive control, but after exact single-site subtraction "
            "it does not robustly change the branch nonlinear interaction law."
        )
    else:
        classification = "METRIC_AUDIT_PASSIVE_POSITIVE_CONTROL_FAILED"
        interpretation = (
            "The locked section-length edit did not produce the required passive "
            "site change, so the nonlinear comparison is not interpretable."
        )

    summary = {
        "object": (
            "direct causal link from intrinsic branch cable length to released-"
            "model nonlinear compact-branch interaction"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "same_gate20_six_branch_basis": True,
            "length_scale": LENGTH_SCALE,
            "diameter_changed": False,
            "topology_changed": False,
            "normalized_synapse_x_changed": False,
            "nmda_kinetics_changed": False,
            "multiplicity_per_site": MULTIPLICITY,
            "thresholds_locked_before_run": {
                "median_nonlinear_effect_min": "log(1.05)",
                "fraction_branches_nonlinear_effect_over_5pct_min": 4.0 / 6.0,
                "median_passive_effect_min": "log(1.01)",
                "spike_guard_fraction_max": 0.0,
            },
        },
        "aggregate": aggregate,
        "branches": rows_out,
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not tune stretch magnitude or synaptic multiplicity after this "
            "result. The audit either connects intrinsic metric to the nonlinear "
            "subunit at the locked 20% intervention or it does not."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print()
    print("Operaattori metric -> nonlinear branch audit")
    print()
    print(
        "median passive site factor:          "
        f"{aggregate['median_passive_site_effect_factor']:.4f}x"
    )
    print(
        "passive branches >1%:                "
        f"{aggregate['fraction_passive_effect_over_1pct']:.3f}"
    )
    print(
        "median nonlinear interaction factor: "
        f"{aggregate['median_nonlinear_interaction_effect_factor']:.4f}x"
    )
    print(
        "nonlinear branches >5%:              "
        f"{aggregate['fraction_branches_nonlinear_effect_over_5pct']:.3f}"
    )
    print(
        "spike guard fraction:                "
        f"{aggregate['fraction_runs_spike_guard']:.3f}"
    )
    print(f"classification: {classification}")

    assert np.all(np.isfinite(nonlinear))
    assert np.all(np.isfinite(passive))


if __name__ == "__main__":
    main()
