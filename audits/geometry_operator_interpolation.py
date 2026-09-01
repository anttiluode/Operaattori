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
    MULTIPLICITY,
    check_branch_identity,
    make_cell,
    recover_branches,
    trace_metrics,
)
from green_circuit import (
    LinearTransport,
    assert_same_kinetics,
    impulse_green_matrix,
    matrix_nrmse,
    measure_raw_conductance,
    noinput_local_trace,
)
from temporal_green_circuit import (
    TIMING_PROGRAMS,
    reduced_solve_timed,
    timed_cluster_trace,
    timed_conductances,
)


ENDPOINTS = (0.80, 1.20)
TARGETS = (0.90, 1.00, 1.10)


def measure_pack(
    fci_root: Path,
    branch: dict,
    kinetics: dict,
) -> tuple[dict, object, object]:
    cell, syn = make_cell(fci_root)
    rows = dendritic_rows(syn)
    sites = np.asarray(branch["sites"], dtype=int)
    check_branch_identity(
        syn, sites, branch["canonical_section"]
    )
    assert_same_kinetics(
        syn, sites, kinetics["tau"]
    )
    return {
        "cell": cell,
        "syn": syn,
        "rows": rows,
        "sites": sites,
    }, cell, syn


def apply_scale(
    syn,
    sites: np.ndarray,
    original_length_um: float,
    scale: float,
) -> float:
    sec = syn.iloc[int(sites[0])]["segments"].sec
    old = float(sec.L)
    if abs(old - float(original_length_um)) > 1e-6:
        raise RuntimeError("fresh model section length differs")
    sec.L = old * float(scale)
    return float(sec.L)


def endpoint_pack(
    fci_root: Path,
    branch: dict,
    scale: float,
    kinetics: dict,
) -> dict:
    state, cell, syn = measure_pack(
        fci_root, branch, kinetics
    )
    sites = state["sites"]
    new_length = apply_scale(
        syn,
        sites,
        branch["section_length_um"],
        scale,
    )
    settle_baselines(syn, state["rows"])
    noinput = noinput_local_trace(
        cell, syn, sites
    )
    local_h, soma_h = impulse_green_matrix(
        cell, syn, sites, noinput
    )
    return {
        "scale": float(scale),
        "section_length_um": new_length,
        "baseline_local_mV": noinput["local"].copy(),
        "local_h": local_h.copy(),
        "soma_h": soma_h.copy(),
        "time": noinput["t"].copy(),
    }


def interpolate_pack(
    lower: dict,
    upper: dict,
    target_scale: float,
) -> dict:
    lo = float(lower["scale"])
    hi = float(upper["scale"])
    alpha = (float(target_scale) - lo) / (hi - lo)
    if not (0.0 <= alpha <= 1.0):
        raise ValueError("target lies outside endpoint interval")
    return {
        "scale": float(target_scale),
        "alpha": float(alpha),
        "baseline_local_mV": (
            (1.0 - alpha) * lower["baseline_local_mV"]
            + alpha * upper["baseline_local_mV"]
        ),
        "local_h": (
            (1.0 - alpha) * lower["local_h"]
            + alpha * upper["local_h"]
        ),
        "soma_h": (
            (1.0 - alpha) * lower["soma_h"]
            + alpha * upper["soma_h"]
        ),
        "time": lower["time"],
    }


def nearest_endpoint(
    endpoint_packs: dict[float, dict],
    target_scale: float,
) -> dict:
    # Locked lower-end tie break for lambda*=1.00.
    if float(target_scale) <= 1.00:
        return endpoint_packs[0.80]
    return endpoint_packs[1.20]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/operator_factorization/"
            "geometry_operator_interpolation.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source is not pinned")

    ref_cell, ref_syn = make_cell(fci_root)
    ref_rows, branches = recover_branches(
        ref_cell, ref_syn, args.cluster_span_um
    )
    settle_baselines(ref_syn, ref_rows)
    ref_site = int(branches[0]["sites"][1])
    kinetics = measure_raw_conductance(
        ref_cell, ref_syn, ref_site
    )
    g_ampa_ref = np.asarray(
        kinetics["g_ampa_uS"], dtype=float
    )
    g_nmda_ref = np.asarray(
        kinetics["g_nmda_raw_uS"], dtype=float
    )

    cases = []
    endpoint_receipts = {}
    interp_converged = []
    direct_converged = []
    spike_flags = []

    for bi, branch in enumerate(branches):
        section = branch["canonical_section"]
        endpoint_packs = {
            scale: endpoint_pack(
                fci_root,
                branch,
                scale,
                kinetics,
            )
            for scale in ENDPOINTS
        }
        endpoint_receipts[section] = {
            str(scale): {
                "section_length_um": endpoint_packs[
                    scale
                ]["section_length_um"]
            }
            for scale in ENDPOINTS
        }

        lower = endpoint_packs[0.80]
        upper = endpoint_packs[1.20]
        if not np.allclose(
            lower["time"],
            upper["time"],
            rtol=0,
            atol=1e-12,
        ):
            raise RuntimeError(
                "endpoint operator time grids differ"
            )

        for target_scale in TARGETS:
            predicted_pack = interpolate_pack(
                lower,
                upper,
                target_scale,
            )
            pred_transport = LinearTransport(
                predicted_pack["local_h"],
                predicted_pack["soma_h"],
            )

            nearest = nearest_endpoint(
                endpoint_packs,
                target_scale,
            )
            nearest_transport = LinearTransport(
                nearest["local_h"],
                nearest["soma_h"],
            )

            # Target full model and directly measured target operators are
            # diagnostics only. They are constructed after the interpolated
            # operator has been fixed from the two endpoints.
            state, cell, syn = measure_pack(
                fci_root,
                branch,
                kinetics,
            )
            sites = state["sites"]
            target_length = apply_scale(
                syn,
                sites,
                branch["section_length_um"],
                target_scale,
            )
            settle_baselines(syn, state["rows"])
            target_noinput = noinput_local_trace(
                cell, syn, sites
            )
            target_local_h, target_soma_h = (
                impulse_green_matrix(
                    cell,
                    syn,
                    sites,
                    target_noinput,
                )
            )
            target_transport = LinearTransport(
                target_local_h,
                target_soma_h,
            )

            if not np.allclose(
                predicted_pack["time"],
                target_noinput["t"],
                rtol=0,
                atol=1e-12,
            ):
                raise RuntimeError(
                    "interpolated/target time grids differ"
                )

            for timing_name, delays in TIMING_PROGRAMS.items():
                actual = timed_cluster_trace(
                    cell,
                    syn,
                    sites,
                    delays,
                    {
                        "t": target_noinput["t"],
                        "soma": target_noinput["soma"],
                    },
                )
                ga, gn = timed_conductances(
                    g_ampa_ref,
                    g_nmda_ref,
                    delays,
                )

                interpolated = reduced_solve_timed(
                    predicted_pack["baseline_local_mV"],
                    pred_transport,
                    ga,
                    gn,
                )
                nearest_result = reduced_solve_timed(
                    nearest["baseline_local_mV"],
                    nearest_transport,
                    ga,
                    gn,
                )
                direct = reduced_solve_timed(
                    target_noinput["local"],
                    target_transport,
                    ga,
                    gn,
                )

                transport_oracle_soma = (
                    target_transport.soma(
                        actual["site_inward_current_nA"]
                    )
                )

                interp_metrics = trace_metrics(
                    actual["soma_depol"],
                    interpolated["soma_depol_mV"],
                    actual["t"],
                )
                nearest_metrics = trace_metrics(
                    actual["soma_depol"],
                    nearest_result["soma_depol_mV"],
                    actual["t"],
                )
                direct_metrics = trace_metrics(
                    actual["soma_depol"],
                    direct["soma_depol_mV"],
                    actual["t"],
                )
                transport_oracle_metrics = trace_metrics(
                    actual["soma_depol"],
                    transport_oracle_soma,
                    actual["t"],
                )
                current_err = matrix_nrmse(
                    actual["site_inward_current_nA"],
                    interpolated["current_nA"],
                )

                interp_converged.append(
                    interpolated["converged"]
                )
                direct_converged.append(
                    direct["converged"]
                )
                spike_flags.append(
                    actual["spike_guard"]
                )

                cases.append(
                    {
                        "branch_index": int(bi),
                        "section": section,
                        "target_scale": float(
                            target_scale
                        ),
                        "target_section_length_um": float(
                            target_length
                        ),
                        "interpolation_alpha": float(
                            predicted_pack["alpha"]
                        ),
                        "nearest_endpoint_scale": float(
                            nearest["scale"]
                        ),
                        "timing_program": timing_name,
                        "delays_ms": list(delays),
                        "interpolated_reduced_soma": interp_metrics,
                        "nearest_endpoint_reduced_soma": nearest_metrics,
                        "direct_target_reduced_oracle": direct_metrics,
                        "target_transport_oracle": (
                            transport_oracle_metrics
                        ),
                        "interpolated_current_nrmse": float(
                            current_err
                        ),
                        "interpolated_fixed_point_converged": bool(
                            interpolated["converged"]
                        ),
                        "direct_fixed_point_converged": bool(
                            direct["converged"]
                        ),
                        "spike_guard": bool(
                            actual["spike_guard"]
                        ),
                    }
                )

                print(
                    f"target={target_scale:.2f} [{bi+1}/6] "
                    f"{section} {timing_name} "
                    f"interp={interp_metrics['nrmse']:.4f} "
                    f"nearest={nearest_metrics['nrmse']:.4f} "
                    f"direct={direct_metrics['nrmse']:.4f} "
                    f"oracle={transport_oracle_metrics['nrmse']:.4f} "
                    f"J={current_err:.4f}"
                )

    interp_err = np.asarray(
        [
            x["interpolated_reduced_soma"]["nrmse"]
            for x in cases
        ],
        dtype=float,
    )
    nearest_err = np.asarray(
        [
            x["nearest_endpoint_reduced_soma"]["nrmse"]
            for x in cases
        ],
        dtype=float,
    )
    direct_err = np.asarray(
        [
            x["direct_target_reduced_oracle"]["nrmse"]
            for x in cases
        ],
        dtype=float,
    )
    transport_oracle_err = np.asarray(
        [
            x["target_transport_oracle"]["nrmse"]
            for x in cases
        ],
        dtype=float,
    )
    current_err = np.asarray(
        [
            x["interpolated_current_nrmse"]
            for x in cases
        ],
        dtype=float,
    )

    scale_medians = {}
    for scale in TARGETS:
        vals = [
            x["interpolated_reduced_soma"]["nrmse"]
            for x in cases
            if abs(x["target_scale"] - scale) < 1e-12
        ]
        scale_medians[str(scale)] = float(
            np.median(vals)
        )

    timing_medians = {}
    for timing_name in TIMING_PROGRAMS:
        vals = [
            x["interpolated_reduced_soma"]["nrmse"]
            for x in cases
            if x["timing_program"] == timing_name
        ]
        timing_medians[timing_name] = float(
            np.median(vals)
        )

    aggregate = {
        "branches": BRANCHES,
        "endpoint_scales": list(ENDPOINTS),
        "target_scales": list(TARGETS),
        "timing_programs": list(
            TIMING_PROGRAMS.keys()
        ),
        "cases": int(len(cases)),
        "multiplicity_per_site": MULTIPLICITY,
        "median_target_transport_oracle_soma_nrmse": float(
            np.median(transport_oracle_err)
        ),
        "median_direct_target_reduced_oracle_soma_nrmse": float(
            np.median(direct_err)
        ),
        "median_interpolated_reduced_soma_nrmse": float(
            np.median(interp_err)
        ),
        "target_scale_median_interpolated_soma_nrmse": (
            scale_medians
        ),
        "timing_median_interpolated_soma_nrmse": (
            timing_medians
        ),
        "median_interpolated_current_nrmse": float(
            np.median(current_err)
        ),
        "median_nearest_endpoint_soma_nrmse": float(
            np.median(nearest_err)
        ),
        "interpolated_to_nearest_median_error_ratio": float(
            np.median(interp_err)
            / (np.median(nearest_err) + 1e-30)
        ),
        "fraction_interpolation_beats_nearest": float(
            np.mean(interp_err < nearest_err)
        ),
        "fraction_interpolated_fixed_points_converged": float(
            np.mean(interp_converged)
        ),
        "fraction_direct_fixed_points_converged": float(
            np.mean(direct_converged)
        ),
        "fraction_spike_guard": float(
            np.mean(spike_flags)
        ),
    }

    target_reduction_ok = (
        aggregate[
            "median_target_transport_oracle_soma_nrmse"
        ] <= 0.01
        and aggregate[
            "median_direct_target_reduced_oracle_soma_nrmse"
        ] <= 0.02
        and aggregate[
            "fraction_direct_fixed_points_converged"
        ] == 1.0
    )

    interpolation_ok = (
        aggregate[
            "median_interpolated_reduced_soma_nrmse"
        ] <= 0.02
        and all(
            v <= 0.03
            for v in scale_medians.values()
        )
        and all(
            v <= 0.03
            for v in timing_medians.values()
        )
        and aggregate[
            "median_interpolated_current_nrmse"
        ] <= 0.02
        and aggregate[
            "interpolated_to_nearest_median_error_ratio"
        ] <= 0.75
        and aggregate[
            "fraction_interpolation_beats_nearest"
        ] >= 48.0 / 72.0
        and aggregate[
            "fraction_interpolated_fixed_points_converged"
        ] == 1.0
        and aggregate["fraction_spike_guard"] == 0.0
    )

    if target_reduction_ok and interpolation_ok:
        classification = (
            "GEOMETRY_INTERPOLATES_REDUCED_OPERATOR"
        )
        interpretation = (
            "Linear interpolation of the complete local Green matrix, soma "
            "transport bank and no-input local baseline between 0.8x and 1.2x "
            "branch-length endpoints predicts unseen intermediate geometries "
            "across asynchronous inputs without target-geometry operator "
            "measurement."
        )
    elif target_reduction_ok:
        classification = (
            "REDUCED_CIRCUIT_VALID_OPERATOR_GEOMETRY_MAP_NOT_LINEAR"
        )
        interpretation = (
            "Directly measured target operators still reduce the full model, "
            "but simple endpoint interpolation is not accurate enough to map "
            "geometry to operator under the locked rulers."
        )
    else:
        classification = (
            "TARGET_GEOMETRY_REDUCTION_FAILED"
        )
        interpretation = (
            "The directly measured reduced operators fail on the new target "
            "geometry panel, so interpolation cannot be interpreted."
        )

    summary = {
        "object": (
            "held-out geometry-to-operator interpolation for the reduced "
            "nonlinear Green circuit"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "endpoint_scales": list(ENDPOINTS),
            "target_scales": list(TARGETS),
            "interpolation": (
                "sample-wise linear interpolation of no-input local baseline, "
                "3x3 local Green kernels and site-to-soma kernels"
            ),
            "timing_programs": {
                k: list(v)
                for k, v in TIMING_PROGRAMS.items()
            },
            "nearest_endpoint_tie_break_at_1p0": 0.80,
            "no_target_operator_used_for_prediction": True,
            "thresholds_locked_before_run": {
                "target_transport_oracle_median_nrmse_max": 0.01,
                "direct_target_reduced_median_nrmse_max": 0.02,
                "interpolated_median_nrmse_max": 0.02,
                "each_target_scale_interpolated_median_nrmse_max": 0.03,
                "each_timing_interpolated_median_nrmse_max": 0.03,
                "interpolated_current_median_nrmse_max": 0.02,
                "interpolated_to_nearest_median_ratio_max": 0.75,
                "interpolation_win_fraction_min": 48.0 / 72.0,
                "interpolated_convergence_fraction_min": 1.0,
                "spike_guard_fraction_max": 0.0,
            },
        },
        "endpoint_receipts": endpoint_receipts,
        "aggregate": aggregate,
        "cases": cases,
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not add polynomial interpolation, move endpoints, fit target "
            "operators, align traces, weight branches or change target scales "
            "after seeing this result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori geometry -> operator interpolation audit")
    print()
    print(
        "target transport oracle NRMSE:      "
        f"{aggregate['median_target_transport_oracle_soma_nrmse']:.4f}"
    )
    print(
        "direct target reduced NRMSE:        "
        f"{aggregate['median_direct_target_reduced_oracle_soma_nrmse']:.4f}"
    )
    print(
        "interpolated reduced NRMSE:         "
        f"{aggregate['median_interpolated_reduced_soma_nrmse']:.4f}"
    )
    print(
        "interpolated current NRMSE:         "
        f"{aggregate['median_interpolated_current_nrmse']:.4f}"
    )
    print(
        "nearest endpoint NRMSE:             "
        f"{aggregate['median_nearest_endpoint_soma_nrmse']:.4f}"
    )
    print(
        "interpolated/nearest:               "
        f"{aggregate['interpolated_to_nearest_median_error_ratio']:.4f}"
    )
    print(
        "interpolation beats nearest:        "
        f"{aggregate['fraction_interpolation_beats_nearest']:.3f}"
    )
    for scale, value in scale_medians.items():
        print(
            f"target scale {scale:>4s} median:       "
            f"{value:.4f}"
        )
    for name, value in timing_medians.items():
        print(
            f"timing {name:12s} median:       "
            f"{value:.4f}"
        )
    print(
        "interpolated convergence:           "
        f"{aggregate['fraction_interpolated_fixed_points_converged']:.3f}"
    )
    print(
        "spike guard:                        "
        f"{aggregate['fraction_spike_guard']:.3f}"
    )
    print(f"classification: {classification}")

    assert len(cases) == 72
    assert np.all(np.isfinite(interp_err))
    assert np.all(np.isfinite(nearest_err))
    assert np.all(np.isfinite(direct_err))
    assert np.all(np.isfinite(transport_oracle_err))


if __name__ == "__main__":
    main()
