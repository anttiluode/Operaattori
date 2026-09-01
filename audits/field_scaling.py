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
    section_groups,
)
from gate24_real_field_nmda import (
    BEND_DEG,
    FIELD_V_PER_M,
    MULTIPLICITY,
    all_segment_records,
    apical_axis,
    bent_coordinates,
    field_coefficients,
    interaction_for,
    load_fielddrive,
    section_descendants,
    settle_baselines,
    xyz_at,
)


AMPLITUDES_V_PER_M = np.asarray([0.25, 0.50, 1.00, 2.00], dtype=float)
FREQUENCY_HZ = 15.0
R2_MIN = 0.98
EFFECT_TARGET_LOG = math.log(1.05)


def recover_proximal_arm(syn_df, branches: int, sites: int, cluster_span_um: float):
    rows = dendritic_rows(syn_df)
    groups = section_groups(syn_df, rows)
    candidates = []

    for name, row_ids in groups.items():
        clustered, span = compact_midpoint_rows(
            syn_df, row_ids, sites, cluster_span_um
        )
        if len(clustered) != sites:
            continue
        sec = syn_df.iloc[int(clustered[0])]["segments"].sec
        desc = section_descendants(sec)
        candidates.append(
            {
                "section": sec,
                "section_name": name,
                "section_length_um": float(sec.L),
                "sites": np.asarray(clustered, dtype=int),
                "span_um": float(span),
                "descendant_sections": int(len(desc)),
                "descendant_cable_um": float(sum(float(s.L) for s in desc)),
            }
        )

    candidates.sort(reverse=True, key=lambda x: x["section_length_um"])
    candidates = candidates[:branches]
    if len(candidates) != branches:
        raise RuntimeError("could not recover Gate-24 six-branch basis")

    return rows, max(candidates, key=lambda x: x["descendant_cable_um"])


def through_origin_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    denom = float(np.dot(x, x))
    if denom <= 0:
        raise ValueError("zero x norm")
    slope = float(np.dot(x, y) / denom)
    pred = slope * x
    sse = float(np.sum((y - pred) ** 2))
    sst0 = float(np.sum(y ** 2))
    r2 = 1.0 - sse / (sst0 + 1e-30)
    return slope, r2


def sign_consistent(y: np.ndarray, tol: float = 1e-12) -> bool:
    nz = y[np.abs(y) > tol]
    if len(nz) == 0:
        return False
    return bool(np.all(nz > 0) or np.all(nz < 0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=6)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/field_scaling/field_scaling_audit.json"),
    )
    args = ap.parse_args()

    if args.branches != 6 or args.sites != 3:
        raise ValueError("field-scaling audit is locked to Gate-24 six x three basis")
    if git_head(args.fci_root.resolve()) != FCI_COMMIT:
        raise RuntimeError("FCI source is not pinned")

    builder = import_builder(args.fci_root.resolve())
    cell, syn_df = builder.create_cell(
        path=str(args.fci_root.resolve() / MODEL_REL) + "/"
    )
    load_fielddrive(cell)

    rows, arm = recover_proximal_arm(
        syn_df, args.branches, args.sites, args.cluster_span_um
    )
    records = all_segment_records(cell)
    soma_xyz = xyz_at(cell.soma[0], 0.5)
    field_axis = apical_axis(records, soma_xyz)
    original_coords = np.stack([r["xyz"] for r in records], axis=0)
    bent_coords, bend_meta = bent_coordinates(
        records, arm["section"], field_axis, BEND_DEG
    )

    # gate24 field_coefficients is locked to FIELD_V_PER_M == 1.0. Build the
    # unit-field spatial coefficients once and scale only their amplitude.
    coeff_original_unit = field_coefficients(
        original_coords, soma_xyz, field_axis
    ) / FIELD_V_PER_M
    coeff_bent_unit = field_coefficients(
        bent_coords, soma_xyz, field_axis
    ) / FIELD_V_PER_M

    baseline_by_row = settle_baselines(syn_df, rows)
    sites = arm["sites"]

    samples = []
    signed = []
    spike_flags = []

    for amp in AMPLITUDES_V_PER_M:
        coeff_original = coeff_original_unit * float(amp)
        coeff_bent = coeff_bent_unit * float(amp)

        original = interaction_for(
            cell,
            syn_df,
            sites,
            "human",
            baseline_by_row,
            records,
            coeff_original,
            FREQUENCY_HZ,
            field_enabled=True,
        )
        bent = interaction_for(
            cell,
            syn_df,
            sites,
            "human",
            baseline_by_row,
            records,
            coeff_bent,
            FREQUENCY_HZ,
            field_enabled=True,
        )

        io = float(original["local_interaction_ratio"])
        ib = float(bent["local_interaction_ratio"])
        effect = float(math.log((ib + 1e-30) / (io + 1e-30)))
        signed.append(effect)
        spike_flags.extend([original["spike_guard"], bent["spike_guard"]])

        samples.append(
            {
                "field_V_per_m": float(amp),
                "original_interaction_ratio": io,
                "bent_interaction_ratio": ib,
                "signed_log_effect": effect,
                "abs_log_effect": float(abs(effect)),
                "effect_factor": float(math.exp(abs(effect))),
                "original_field_only_local_peak_to_peak_mV": float(
                    original["field_only_local_peak_to_peak_mV"]
                ),
                "bent_field_only_local_peak_to_peak_mV": float(
                    bent["field_only_local_peak_to_peak_mV"]
                ),
                "spike_guard": bool(
                    original["spike_guard"] or bent["spike_guard"]
                ),
            }
        )

    y = np.asarray(signed, dtype=float)
    slope, r2 = through_origin_fit(AMPLITUDES_V_PER_M, y)
    same_sign = sign_consistent(y)
    spike_fraction = float(np.mean(spike_flags))
    linear_ok = bool(r2 >= R2_MIN and same_sign and spike_fraction == 0.0)

    if linear_ok and abs(slope) > 1e-15:
        field_for_5pct = float(EFFECT_TARGET_LOG / abs(slope))
    else:
        field_for_5pct = None

    if linear_ok:
        classification = "WEAK_FIELD_BEND_EFFECT_IN_LOCAL_LINEAR_REGIME"
        interpretation = (
            "Across the locked 0.25-2 V/m panel, the signed geometry-dependent "
            "nonlinear interaction effect is consistent with a through-origin "
            "small-signal response. The 5% field is only a local first-order "
            "scale estimate, not a prediction that the response remains linear "
            "at that amplitude."
        )
    else:
        classification = "NO_RELIABLE_LOCAL_LINEAR_FIELD_SCALING"
        interpretation = (
            "The locked weak-field panel does not support a stable through-origin "
            "linear extrapolation of the geometry-dependent nonlinear effect."
        )

    summary = {
        "object": (
            "bounded small-signal amplitude audit of Gate-24 proximal "
            "real-field bend effect"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "branch": "proximal_large_subtree",
            "section": arm["section_name"],
            "descendant_cable_um": arm["descendant_cable_um"],
            "sites": sites.tolist(),
            "compact_span_um": arm["span_um"],
            "bend_degrees": BEND_DEG,
            "frequency_hz": FREQUENCY_HZ,
            "field_amplitudes_V_per_m": AMPLITUDES_V_PER_M.tolist(),
            "multiplicity_per_site": MULTIPLICITY,
            "virtual_synapses": int(len(sites) * MULTIPLICITY),
            "only_independent_variable": "field amplitude",
            "fit": "signed log bend effect through physical origin",
            "r2_min_for_extrapolation": R2_MIN,
            "target_effect_log": EFFECT_TARGET_LOG,
        },
        "bend": bend_meta,
        "samples": samples,
        "aggregate": {
            "signed_effects": y.tolist(),
            "through_origin_slope_log_per_V_per_m": slope,
            "through_origin_r2": r2,
            "sign_consistent": same_sign,
            "fraction_runs_spike_guard": spike_fraction,
            "local_extrapolated_field_for_5pct_V_per_m": field_for_5pct,
            "extrapolation_is_local_only": True,
        },
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not add larger field amplitudes after seeing this audit. If the "
            "local estimate lies beyond the audited panel, report that scale "
            "rather than simulating upward until the 5% ruler is crossed."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("Operaattori field-scaling audit")
    print()
    print(f"section:                         {arm['section_name']}")
    print(f"max bend displacement:           {bend_meta['max_segment_displacement_um']:.3f} um")
    for row in samples:
        print(
            f"E={row['field_V_per_m']:.2f} V/m  "
            f"signed_log={row['signed_log_effect']:+.6e}  "
            f"factor={row['effect_factor']:.6f}x"
        )
    print()
    print(f"through-origin slope:            {slope:+.6e} log/(V/m)")
    print(f"through-origin R2:               {r2:.6f}")
    print(f"sign consistent:                 {same_sign}")
    print(f"spike guard fraction:            {spike_fraction:.3f}")
    if field_for_5pct is None:
        print("local field estimate for 5%:     unavailable")
    else:
        print(f"local field estimate for 5%:     {field_for_5pct:.3f} V/m")
    print(f"classification:                   {classification}")

    assert np.all(np.isfinite(y))
    assert spike_fraction == 0.0


if __name__ == "__main__":
    main()
