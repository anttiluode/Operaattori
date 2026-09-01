from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

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
    section_groups,
)
from gate17_superposition_attack import block


RAT_GAMMA = 0.062
HUMAN_RATIO = 0.00131 / 0.00088
FIELD_V_PER_M = 1.0
FREQUENCIES_HZ = (10.0, 15.0, 20.0)
BEND_DEG = 35.0
MULTIPLICITY = 8.0
EVENT_MS = 300.0
TSTOP_MS = 410.0
DT_MS = 0.025


def xyz_at(sec, x: float) -> np.ndarray:
    from neuron import h

    n = int(h.n3d(sec=sec))
    if n < 2:
        raise RuntimeError(f"section {sec.name()} lacks pt3d geometry")
    arcs = np.asarray([float(h.arc3d(i, sec=sec)) for i in range(n)])
    xyz = np.asarray(
        [
            [
                float(h.x3d(i, sec=sec)),
                float(h.y3d(i, sec=sec)),
                float(h.z3d(i, sec=sec)),
            ]
            for i in range(n)
        ],
        dtype=float,
    )
    target = float(x) * float(sec.L)
    if target <= arcs[0]:
        return xyz[0].copy()
    if target >= arcs[-1]:
        return xyz[-1].copy()
    j = int(np.searchsorted(arcs, target, side="right"))
    a0, a1 = float(arcs[j - 1]), float(arcs[j])
    w = (target - a0) / (a1 - a0 + 1e-30)
    return (1.0 - w) * xyz[j - 1] + w * xyz[j]


def section_descendants(sec) -> list:
    from neuron import h

    out = []
    stack = [sec]
    seen = set()
    while stack:
        s = stack.pop()
        name = s.name()
        if name in seen:
            continue
        seen.add(name)
        out.append(s)
        ref = h.SectionRef(sec=s)
        for child in ref.child:
            stack.append(child)
    return out


def all_segment_records(cell):
    out = []
    for sec in cell.all:
        for seg in sec:
            out.append(
                {
                    "sec": sec,
                    "sec_name": sec.name(),
                    "x": float(seg.x),
                    "seg": seg,
                    "xyz": xyz_at(sec, float(seg.x)),
                }
            )
    return out


def apical_axis(records, soma_xyz: np.ndarray) -> np.ndarray:
    coords = np.stack(
        [
            r["xyz"]
            for r in records
            if "apic" in r["sec_name"]
        ],
        axis=0,
    )
    X = coords - soma_xyz[None, :]
    cov = X.T @ X
    vals, vecs = np.linalg.eigh(cov)
    axis = vecs[:, int(np.argmax(vals))]
    mean = np.mean(X, axis=0)
    if float(np.dot(axis, mean)) < 0:
        axis = -axis
    axis = axis / np.linalg.norm(axis)
    return axis


def safe_rotation_axis(field_axis: np.ndarray, subtree_axis: np.ndarray) -> np.ndarray:
    axis = np.cross(subtree_axis, field_axis)
    n = float(np.linalg.norm(axis))
    if n < 1e-9:
        trial = np.asarray([1.0, 0.0, 0.0])
        if abs(float(np.dot(trial, field_axis))) > 0.9:
            trial = np.asarray([0.0, 1.0, 0.0])
        axis = np.cross(field_axis, trial)
        n = float(np.linalg.norm(axis))
    return axis / n


def rodrigues(axis: np.ndarray, angle_deg: float) -> np.ndarray:
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    x, y, z = a
    K = np.asarray(
        [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]],
        dtype=float,
    )
    th = math.radians(float(angle_deg))
    I = np.eye(3)
    return I + math.sin(th) * K + (1.0 - math.cos(th)) * (K @ K)


def bent_coordinates(
    records,
    bend_sec,
    field_axis: np.ndarray,
    angle_deg: float,
) -> tuple[np.ndarray, dict]:
    descendants = section_descendants(bend_sec)
    names = {s.name() for s in descendants}
    pivot = xyz_at(bend_sec, 0.0)

    moved_xyz = np.stack(
        [r["xyz"] for r in records if r["sec_name"] in names],
        axis=0,
    )
    if len(moved_xyz) == 0:
        raise RuntimeError("empty bend subtree")
    subtree_axis = np.mean(moved_xyz, axis=0) - pivot
    if np.linalg.norm(subtree_axis) < 1e-9:
        subtree_axis = xyz_at(bend_sec, 1.0) - pivot

    rot_axis = safe_rotation_axis(field_axis, subtree_axis)
    R = rodrigues(rot_axis, angle_deg)

    coords = np.stack([r["xyz"] for r in records], axis=0)
    mask = np.asarray(
        [r["sec_name"] in names for r in records],
        dtype=bool,
    )
    rel = coords[mask] - pivot[None, :]
    coords[mask] = pivot[None, :] + rel @ R.T

    displacement = np.linalg.norm(
        coords[mask] - moved_xyz,
        axis=1,
    )
    meta = {
        "section": bend_sec.name(),
        "descendant_sections": int(len(descendants)),
        "descendant_cable_um": float(sum(float(s.L) for s in descendants)),
        "moved_segments": int(np.sum(mask)),
        "max_segment_displacement_um": float(np.max(displacement)),
        "mean_segment_displacement_um": float(np.mean(displacement)),
        "rotation_axis": rot_axis.tolist(),
        "pivot_xyz_um": pivot.tolist(),
    }
    return coords, meta


def field_coefficients(
    coords: np.ndarray,
    soma_xyz: np.ndarray,
    field_axis: np.ndarray,
) -> np.ndarray:
    # Vext[mV] = -E[V/m] * projection[um] * 1e-3.
    proj_um = (coords - soma_xyz[None, :]) @ field_axis
    return -1e-3 * FIELD_V_PER_M * proj_um


def load_fielddrive(cell):
    import neuron
    from neuron import h

    # CI compiles mechanisms/ into the repo root x86_64 directory.
    neuron.load_mechanisms(str(ROOT))
    for sec in cell.all:
        sec.insert("extracellular")
        sec.insert("fielddrive")
    for sec in cell.all:
        for seg in sec:
            h.setpointer(
                seg._ref_e_extracellular,
                "ex",
                seg.fielddrive,
            )


def set_coefficients(records, coeff: np.ndarray) -> None:
    if len(records) != len(coeff):
        raise ValueError("coefficient length mismatch")
    for r, c in zip(records, coeff):
        r["seg"].fielddrive.coeff = float(c)


def field_waveform(freq_hz: float, *, enabled: bool = True):
    from neuron import h

    t = np.arange(0.0, TSTOP_MS + DT_MS * 0.5, DT_MS)
    if enabled:
        # Same external field for ORIGINAL and BENT. Positive peak at EVENT_MS.
        drive = np.cos(
            2.0 * np.pi * float(freq_hz) * (t - EVENT_MS) / 1000.0
        )
    else:
        drive = np.zeros_like(t)
    tv = h.Vector(t)
    dv = h.Vector(drive)
    dv.play(h._ref_drive_fielddrive, tv, 1)
    return tv, dv


def settle_baselines(syn_df, rows: np.ndarray) -> dict[int, float]:
    from neuron import h

    h.drive_fielddrive = 0.0
    h.dt = DT_MS
    h.finitialize(V_INIT_MV)
    h.continuerun(55.0)
    return {
        int(i): float(syn_df.iloc[int(i)]["segments"].v)
        for i in rows
    }


def configure_rows(
    syn_df,
    rows: np.ndarray,
    condition: str,
    baseline_by_row: dict[int, float],
) -> None:
    for i in rows:
        row = syn_df.iloc[int(i)]
        syn = row["exc_synapses"]
        nc = row["exc_netcons"]
        vrest = float(baseline_by_row[int(i)])

        if condition == "human":
            syn.gamma = HUMAN_GAMMA
            syn.NMDA_ratio = HUMAN_RATIO
        elif condition == "gamma062_restmatched":
            syn.gamma = RAT_GAMMA
            syn.NMDA_ratio = (
                HUMAN_RATIO
                * block(vrest, HUMAN_GAMMA)
                / block(vrest, RAT_GAMMA)
            )
        else:
            raise ValueError(condition)

        nc.weight[0] = BASE_WEIGHT_US * MULTIPLICITY


def run_trace(
    cell,
    syn_df,
    target_rows: np.ndarray,
    active_rows: np.ndarray,
    condition: str,
    baseline_by_row: dict[int, float],
    records,
    coeff: np.ndarray,
    freq_hz: float,
    *,
    field_enabled: bool,
) -> dict:
    from neuron import h

    target_rows = np.asarray(target_rows, dtype=int)
    active_rows = np.asarray(active_rows, dtype=int)
    configure_rows(
        syn_df,
        target_rows,
        condition,
        baseline_by_row,
    )
    set_coefficients(records, coeff)
    tv_play, dv_play = field_waveform(freq_hz, enabled=field_enabled)

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
    local_vecs = [
        h.Vector().record(
            syn_df.iloc[int(i)]["segments"]._ref_v
        )
        for i in target_rows
    ]
    nmda_vecs = [
        h.Vector().record(
            syn_df.iloc[int(i)]["exc_synapses"]._ref_i_NMDA
        )
        for i in target_rows
    ]

    h.dt = DT_MS
    h.finitialize(V_INIT_MV)
    h.fcurrent()
    for i in active_rows:
        syn_df.iloc[int(i)]["exc_netcons"].event(EVENT_MS)
    h.continuerun(TSTOP_MS)

    # Keep play vectors alive until after continuerun.
    _ = (tv_play, dv_play)

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    local = np.stack(
        [np.asarray(v, dtype=float) for v in local_vecs],
        axis=0,
    )
    nmda = np.stack(
        [np.asarray(v, dtype=float) for v in nmda_vecs],
        axis=0,
    )
    if not (
        np.all(np.isfinite(soma))
        and np.all(np.isfinite(local))
        and np.all(np.isfinite(nmda))
    ):
        raise FloatingPointError("non-finite Gate-24 trace")

    post = (t >= EVENT_MS) & (t <= EVENT_MS + 90.0)
    pre = (t >= EVENT_MS - 40.0) & (t < EVENT_MS - 5.0)
    return {
        "t": t[post],
        "soma": soma[post],
        "local": local[:, post],
        "nmda": nmda[:, post],
        "pre_soma_peak_mV": float(np.max(soma[pre])),
        "post_soma_peak_mV": float(np.max(soma[post])),
    }


def positive_auc(trace: np.ndarray, t: np.ndarray) -> float:
    return float(np.trapezoid(np.maximum(trace, 0.0), t))


def interaction_for(
    cell,
    syn_df,
    sites: np.ndarray,
    condition: str,
    baseline_by_row: dict[int, float],
    records,
    coeff: np.ndarray,
    freq_hz: float,
    *,
    field_enabled: bool,
) -> dict:
    field_only = run_trace(
        cell,
        syn_df,
        sites,
        np.empty(0, dtype=int),
        condition,
        baseline_by_row,
        records,
        coeff,
        freq_hz,
        field_enabled=field_enabled,
    )
    together = run_trace(
        cell,
        syn_df,
        sites,
        sites,
        condition,
        baseline_by_row,
        records,
        coeff,
        freq_hz,
        field_enabled=field_enabled,
    )
    singles = [
        run_trace(
            cell,
            syn_df,
            sites,
            np.asarray([site], dtype=int),
            condition,
            baseline_by_row,
            records,
            coeff,
            freq_hz,
            field_enabled=field_enabled,
        )
        for site in sites
    ]

    t = together["t"]
    for tr in [field_only, *singles]:
        if not np.allclose(tr["t"], t, rtol=0, atol=1e-12):
            raise RuntimeError("Gate-24 time grids differ")

    actual_local = together["local"] - field_only["local"]
    single_local = [
        tr["local"] - field_only["local"]
        for tr in singles
    ]
    predicted_local = np.sum(
        np.stack(single_local, axis=0),
        axis=0,
    )

    actual_mean = np.mean(actual_local, axis=0)
    predicted_mean = np.mean(predicted_local, axis=0)
    actual_auc = positive_auc(actual_mean, t)
    predicted_auc = positive_auc(predicted_mean, t)
    interaction = float(actual_auc / (predicted_auc + 1e-30))

    actual_soma = together["soma"] - field_only["soma"]
    single_soma = [
        tr["soma"] - field_only["soma"]
        for tr in singles
    ]
    predicted_soma = np.sum(np.stack(single_soma, axis=0), axis=0)
    soma_auc = positive_auc(actual_soma, t)
    soma_pred_auc = positive_auc(predicted_soma, t)

    nmda_charge = float(
        np.trapezoid(
            np.sum(np.abs(together["nmda"]), axis=0),
            t,
        )
    )

    spike_guard = bool(
        together["post_soma_peak_mV"] >= -20.0
        or field_only["post_soma_peak_mV"] >= -20.0
    )

    return {
        "local_interaction_ratio": interaction,
        "actual_local_auc_mV_ms": actual_auc,
        "independent_local_auc_mV_ms": predicted_auc,
        "soma_interaction_ratio": float(
            soma_auc / (soma_pred_auc + 1e-30)
        ),
        "nmda_abs_current_auc_nA_ms": nmda_charge,
        "field_only_local_peak_to_peak_mV": float(
            np.median(
                np.max(field_only["local"], axis=1)
                - np.min(field_only["local"], axis=1)
            )
        ),
        "field_only_soma_peak_to_peak_mV": float(
            np.max(field_only["soma"]) - np.min(field_only["soma"])
        ),
        "spike_guard": spike_guard,
    }


def rel(a: float, b: float) -> float:
    return float(abs(a - b) / (0.5 * (abs(a) + abs(b)) + 1e-30))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--branches", type=int, default=6)
    ap.add_argument("--sites", type=int, default=3)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("results/gate24/gate24_real_field_nmda.json"),
    )
    args = ap.parse_args()

    if args.branches != 6 or args.sites != 3:
        raise ValueError("Gate 24 is locked to the Gate-20 six x three basis")
    if git_head(args.fci_root.resolve()) != FCI_COMMIT:
        raise RuntimeError("FCI source is not pinned")

    builder = import_builder(args.fci_root.resolve())
    cell, syn_df = builder.create_cell(
        path=str(args.fci_root.resolve() / MODEL_REL) + "/"
    )
    load_fielddrive(cell)

    rows = dendritic_rows(syn_df)
    groups = section_groups(syn_df, rows)
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
        desc = section_descendants(sec)
        candidates.append(
            {
                "section": sec,
                "section_name": name,
                "section_length_um": float(sec.L),
                "sites": np.asarray(clustered, dtype=int),
                "span_um": float(span),
                "descendant_sections": int(len(desc)),
                "descendant_cable_um": float(
                    sum(float(s.L) for s in desc)
                ),
            }
        )

    candidates.sort(
        reverse=True,
        key=lambda x: x["section_length_um"],
    )
    candidates = candidates[:6]
    if len(candidates) != 6:
        raise RuntimeError("Gate 24 could not recover Gate-20 branch basis")

    by_extent = sorted(
        candidates,
        key=lambda x: x["descendant_cable_um"],
    )
    arms = {
        "distal_small_subtree": by_extent[0],
        "proximal_large_subtree": by_extent[-1],
    }

    records = all_segment_records(cell)
    soma_xyz = xyz_at(cell.soma[0], 0.5)
    field_axis = apical_axis(records, soma_xyz)
    original_coords = np.stack([r["xyz"] for r in records], axis=0)
    coeff_original = field_coefficients(
        original_coords,
        soma_xyz,
        field_axis,
    )

    baseline_by_row = settle_baselines(syn_df, rows)

    results = {}
    all_spikes = []
    for arm_name, arm in arms.items():
        bent_coords, bend_meta = bent_coordinates(
            records,
            arm["section"],
            field_axis,
            BEND_DEG,
        )
        coeff_bent = field_coefficients(
            bent_coords,
            soma_xyz,
            field_axis,
        )

        sites = arm["sites"]
        arm_out = {
            "section": arm["section_name"],
            "sites": sites.tolist(),
            "compact_span_um": arm["span_um"],
            "section_length_um": arm["section_length_um"],
            "descendant_sections": arm["descendant_sections"],
            "descendant_cable_um": arm["descendant_cable_um"],
            "bend": bend_meta,
            "field_projection": {
                "median_abs_original_coeff_mV": float(
                    np.median(np.abs(coeff_original))
                ),
                "median_abs_bent_minus_original_coeff_mV_all_segments": float(
                    np.median(np.abs(coeff_bent - coeff_original))
                ),
                "max_abs_bent_minus_original_coeff_mV_all_segments": float(
                    np.max(np.abs(coeff_bent - coeff_original))
                ),
                "target_site_coeff_original_mV": [
                    float(
                        -1e-3
                        * FIELD_V_PER_M
                        * np.dot(
                            xyz_at(
                                syn_df.iloc[int(i)]["segments"].sec,
                                float(syn_df.iloc[int(i)]["segments"].x),
                            )
                            - soma_xyz,
                            field_axis,
                        )
                    )
                    for i in sites
                ],
            },
            "frequencies": {},
        }

        # Zero-field geometry control. Since geometry enters only through field
        # coefficients, ORIGINAL and BENT must match when drive=0.
        zero_o = interaction_for(
            cell,
            syn_df,
            sites,
            "human",
            baseline_by_row,
            records,
            coeff_original,
            15.0,
            field_enabled=False,
        )
        zero_b = interaction_for(
            cell,
            syn_df,
            sites,
            "human",
            baseline_by_row,
            records,
            coeff_bent,
            15.0,
            field_enabled=False,
        )
        zero_rel = rel(
            zero_o["local_interaction_ratio"],
            zero_b["local_interaction_ratio"],
        )

        # Material-locked control: "bent" material receives the original
        # coefficients, which is exactly the original physical forcing.
        locked = interaction_for(
            cell,
            syn_df,
            sites,
            "human",
            baseline_by_row,
            records,
            coeff_original,
            15.0,
            field_enabled=True,
        )
        locked_repeat = interaction_for(
            cell,
            syn_df,
            sites,
            "human",
            baseline_by_row,
            records,
            coeff_original,
            15.0,
            field_enabled=True,
        )
        locked_rel = rel(
            locked["local_interaction_ratio"],
            locked_repeat["local_interaction_ratio"],
        )

        arm_out["controls"] = {
            "zero_field_original_vs_bent_relative_difference": zero_rel,
            "material_locked_repeat_relative_difference": locked_rel,
        }

        for freq in FREQUENCIES_HZ:
            frow = {}
            for condition in (
                "human",
                "gamma062_restmatched",
            ):
                original = interaction_for(
                    cell,
                    syn_df,
                    sites,
                    condition,
                    baseline_by_row,
                    records,
                    coeff_original,
                    freq,
                    field_enabled=True,
                )
                bent = interaction_for(
                    cell,
                    syn_df,
                    sites,
                    condition,
                    baseline_by_row,
                    records,
                    coeff_bent,
                    freq,
                    field_enabled=True,
                )
                effect_log = float(
                    abs(
                        math.log(
                            (bent["local_interaction_ratio"] + 1e-30)
                            / (original["local_interaction_ratio"] + 1e-30)
                        )
                    )
                )
                frow[condition] = {
                    "original": original,
                    "bent": bent,
                    "bend_interaction_effect_abs_log": effect_log,
                    "bend_interaction_effect_factor": float(
                        math.exp(effect_log)
                    ),
                }
                all_spikes.extend(
                    [original["spike_guard"], bent["spike_guard"]]
                )
            arm_out["frequencies"][str(freq)] = frow
        results[arm_name] = arm_out
        print(
            f"{arm_name:24s} {arm['section_name']} "
            f"desc={arm['descendant_cable_um']:.1f}um "
            f"moved={bend_meta['max_segment_displacement_um']:.1f}um "
            f"zero={zero_rel:.3e} locked={locked_rel:.3e}"
        )

    def effects(arm_name: str, condition: str) -> np.ndarray:
        return np.asarray(
            [
                results[arm_name]["frequencies"][str(f)][condition][
                    "bend_interaction_effect_abs_log"
                ]
                for f in FREQUENCIES_HZ
            ],
            dtype=float,
        )

    prox_h = effects("proximal_large_subtree", "human")
    dist_h = effects("distal_small_subtree", "human")
    prox_r = effects(
        "proximal_large_subtree",
        "gamma062_restmatched",
    )
    threshold = math.log(1.05)

    max_control = max(
        results[a]["controls"][
            "zero_field_original_vs_bent_relative_difference"
        ]
        for a in results
    )
    max_locked = max(
        results[a]["controls"][
            "material_locked_repeat_relative_difference"
        ]
        for a in results
    )

    aggregate = {
        "field_V_per_m": FIELD_V_PER_M,
        "frequencies_hz": list(FREQUENCIES_HZ),
        "bend_degrees": BEND_DEG,
        "virtual_synapses": int(3 * MULTIPLICITY),
        "max_zero_field_control_relative_difference": float(max_control),
        "max_material_locked_relative_difference": float(max_locked),
        "proximal_human_effect_abs_log": prox_h.tolist(),
        "distal_human_effect_abs_log": dist_h.tolist(),
        "proximal_restmatched_effect_abs_log": prox_r.tolist(),
        "proximal_human_median_effect_abs_log": float(np.median(prox_h)),
        "distal_human_median_effect_abs_log": float(np.median(dist_h)),
        "proximal_human_median_effect_factor": float(
            math.exp(float(np.median(prox_h)))
        ),
        "distal_human_median_effect_factor": float(
            math.exp(float(np.median(dist_h)))
        ),
        "proximal_frequencies_over_5pct": int(np.sum(prox_h >= threshold)),
        "proximal_to_distal_median_effect_ratio": float(
            np.median(prox_h) / (np.median(dist_h) + 1e-30)
        ),
        "fraction_runs_spike_guard": float(np.mean(all_spikes)),
    }

    controls_ok = max_control < 1e-6 and max_locked < 1e-6
    primary_ok = (
        aggregate["proximal_frequencies_over_5pct"] >= 2
        and aggregate["proximal_to_distal_median_effect_ratio"] >= 2.0
        and aggregate["fraction_runs_spike_guard"] == 0.0
    )

    max_field_polarization = max(
        results[a]["frequencies"][str(f)]["human"]["original"][
            "field_only_local_peak_to_peak_mV"
        ]
        for a in results
        for f in FREQUENCIES_HZ
    )

    if controls_ok and primary_ok:
        classification = (
            "REAL_FIELD_BEND_MODULATES_NONLINEAR_COMPARTMENT"
        )
        interpretation = (
            "A rigid re-embedding that leaves the intrinsic cell-1125 cable "
            "unchanged alters the nonlinear compact-branch interaction under "
            "the same physical extracellular field, and the large-subtree "
            "arm exceeds the distal anatomical attacker."
        )
    elif controls_ok and max_field_polarization > 1e-4:
        classification = (
            "FIELD_COUPLING_PRESENT_BUT_NOT_COMPARTMENT_SELECTIVE"
        )
        interpretation = (
            "The physical extracellular field polarizes the released model "
            "and the isometric geometry changes field coefficients, but the "
            "preregistered nonlinear branch-interaction effect does not clear "
            "the large-subtree-versus-distal ruler."
        )
    else:
        classification = "NO_ROBUST_REAL_FIELD_BEND_EFFECT"
        interpretation = (
            "The real-field scaffold assay does not establish a robust "
            "geometry-dependent nonlinear-compartment effect under the "
            "locked 1 V/m, 35-degree, 24-synapse protocol."
        )

    summary = {
        "gate": 24,
        "object": (
            "uniform extracellular field coupled through NEURON to an "
            "isometrically re-embedded real cell-1125 scaffold with released "
            "human NMDA synapses"
        ),
        "fci_commit": FCI_COMMIT,
        "literature_scope": {
            "aspart_2018": "10.1371/journal.pcbi.1006124",
            "fan_2023_2024": "10.1007/s11571-022-09922-y",
            "novelty_claim_field_nmda": False,
            "released_model_dendrites_are_passive": True,
            "active_ih_resonance_claimed": False,
        },
        "protocol": {
            "gate20_branch_basis_reused": True,
            "field_V_per_m": FIELD_V_PER_M,
            "frequencies_hz": list(FREQUENCIES_HZ),
            "bend_degrees": BEND_DEG,
            "multiplicity_per_site": MULTIPLICITY,
            "sites_per_branch": 3,
            "virtual_synapses": int(3 * MULTIPLICITY),
            "event_ms": EVENT_MS,
            "intrinsic_section_geometry_modified": False,
            "extracellular_potential_from_world_coordinates": True,
            "field_axis": field_axis.tolist(),
            "thresholds_locked_before_run": {
                "control_relative_difference_max": 1e-6,
                "proximal_frequencies_over_5pct_min": 2,
                "proximal_to_distal_median_effect_ratio_min": 2.0,
                "spike_guard_fraction_max": 0.0,
            },
        },
        "aggregate": aggregate,
        "arms": results,
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not open Gate 25 in this session. Do not rescue a negative "
            "Gate 24 by increasing field amplitude, scanning bend angle or "
            "tuning synaptic multiplicity."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori Gate 24 — real extracellular field x NMDA compartment")
    print()
    print(
        "zero-field max control:             "
        f"{aggregate['max_zero_field_control_relative_difference']:.3e}"
    )
    print(
        "material-locked max control:         "
        f"{aggregate['max_material_locked_relative_difference']:.3e}"
    )
    print(
        "proximal HUMAN median bend factor:   "
        f"{aggregate['proximal_human_median_effect_factor']:.4f}x"
    )
    print(
        "distal HUMAN median bend factor:     "
        f"{aggregate['distal_human_median_effect_factor']:.4f}x"
    )
    print(
        "proximal frequencies >5%:           "
        f"{aggregate['proximal_frequencies_over_5pct']} / 3"
    )
    print(
        "proximal/distal median effect ratio: "
        f"{aggregate['proximal_to_distal_median_effect_ratio']:.3f}x"
    )
    print(
        "spike guard fraction:               "
        f"{aggregate['fraction_runs_spike_guard']:.3f}"
    )
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert np.all(np.isfinite(prox_h))
    assert np.all(np.isfinite(dist_h))
    assert classification in {
        "REAL_FIELD_BEND_MODULATES_NONLINEAR_COMPARTMENT",
        "FIELD_COUPLING_PRESENT_BUT_NOT_COMPARTMENT_SELECTIVE",
        "NO_ROBUST_REAL_FIELD_BEND_EFFECT",
    }


if __name__ == "__main__":
    main()
