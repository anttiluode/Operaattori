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
    DT_MS,
    EVENT_MS,
    MULTIPLICITY,
    POST_MS,
    TSTOP_MS,
    check_branch_identity,
    configure_human,
    make_cell,
    recover_branches,
    rms,
    trace_metrics,
)
from green_circuit import (
    DAMPING,
    FP_MAX_ITER,
    FP_TOL,
    GEOMETRIES,
    LinearTransport,
    assert_same_kinetics,
    impulse_green_matrix,
    matrix_nrmse,
    measure_raw_conductance,
    noinput_local_trace,
    synapse_law,
)


TIMING_PROGRAMS = {
    "synchronous": (0.0, 0.0, 0.0),
    "forward_5": (0.0, 5.0, 10.0),
    "reverse_5": (10.0, 5.0, 0.0),
    "spread_15": (0.0, 15.0, 30.0),
}


def timed_cluster_trace(
    cell,
    syn_df,
    sites: np.ndarray,
    delays_ms: tuple[float, float, float],
    noinput: dict,
) -> dict:
    from neuron import h

    sites = np.asarray(sites, dtype=int)
    if len(sites) != 3 or len(delays_ms) != 3:
        raise ValueError("temporal audit is locked to three sites")

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
    h.finitialize(-76.0)
    h.fcurrent()
    for i, delay in zip(sites, delays_ms):
        syn_df.iloc[int(i)]["exc_netcons"].event(
            EVENT_MS + float(delay)
        )
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    ampa = np.stack(
        [np.asarray(v, dtype=float) for v in ampa_vecs],
        axis=0,
    )
    nmda = np.stack(
        [np.asarray(v, dtype=float) for v in nmda_vecs],
        axis=0,
    )

    post = (
        (t >= EVENT_MS)
        & (t <= EVENT_MS + POST_MS + 1e-9)
    )
    tp = t[post] - EVENT_MS
    if not np.allclose(
        tp, noinput["t"], rtol=0, atol=1e-12
    ):
        raise RuntimeError("temporal matched-control grid differs")

    soma_dep = soma[post] - noinput["soma"]
    inward = -(ampa[:, post] + nmda[:, post])

    if not (
        np.all(np.isfinite(soma_dep))
        and np.all(np.isfinite(inward))
    ):
        raise FloatingPointError("non-finite temporal full-model trace")

    return {
        "t": tp,
        "soma_depol": soma_dep,
        "site_inward_current_nA": inward,
        "soma_peak_absolute_mV": float(np.max(soma[post])),
        "spike_guard": bool(np.max(soma[post]) >= -20.0),
    }


def shift_template(
    template: np.ndarray,
    delay_ms: float,
) -> np.ndarray:
    template = np.asarray(template, dtype=float)
    steps_float = float(delay_ms) / DT_MS
    steps = int(round(steps_float))
    if abs(steps_float - steps) > 1e-9:
        raise ValueError(
            f"delay {delay_ms} ms is not an integer DT multiple"
        )
    out = np.zeros_like(template)
    if steps < len(template):
        out[steps:] = template[: len(template) - steps]
    return out


def timed_conductances(
    g_ampa_ref: np.ndarray,
    g_nmda_ref: np.ndarray,
    delays_ms: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    ga = np.stack(
        [
            shift_template(g_ampa_ref, delay)
            for delay in delays_ms
        ],
        axis=0,
    )
    gn = np.stack(
        [
            shift_template(g_nmda_ref, delay)
            for delay in delays_ms
        ],
        axis=0,
    )
    return ga, gn


def reduced_solve_timed(
    noinput_local: np.ndarray,
    transport: LinearTransport,
    ga: np.ndarray,
    gn: np.ndarray,
) -> dict:
    noinput_local = np.asarray(noinput_local, dtype=float)
    ga = np.asarray(ga, dtype=float)
    gn = np.asarray(gn, dtype=float)

    current = synapse_law(noinput_local, ga, gn)
    converged = False
    final_err = float("inf")
    iterations = 0

    for iteration in range(1, FP_MAX_ITER + 1):
        voltage = noinput_local + transport.local(current)
        target = synapse_law(voltage, ga, gn)
        updated = (
            (1.0 - DAMPING) * current
            + DAMPING * target
        )
        final_err = rms(updated - current) / (
            rms(updated) + 1e-30
        )
        current = updated
        iterations = iteration
        if final_err <= FP_TOL:
            converged = True
            break

    voltage = noinput_local + transport.local(current)
    soma = transport.soma(current)
    return {
        "current_nA": current,
        "local_voltage_mV": voltage,
        "soma_depol_mV": soma,
        "converged": converged,
        "iterations": int(iterations),
        "final_relative_current_update": float(final_err),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/operator_factorization/temporal_green_circuit.json"
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
    converged_flags = []
    spike_flags = []

    for bi, branch in enumerate(branches):
        section_key = branch["canonical_section"]
        sites_all = np.asarray(branch["sites"], dtype=int)

        for scale in GEOMETRIES:
            cell, syn = make_cell(fci_root)
            rows = dendritic_rows(syn)
            check_branch_identity(
                syn, sites_all, section_key
            )
            assert_same_kinetics(
                syn, sites_all, kinetics["tau"]
            )

            sec = syn.iloc[int(sites_all[0])]["segments"].sec
            old_length = float(sec.L)
            if abs(
                old_length - branch["section_length_um"]
            ) > 1e-6:
                raise RuntimeError("fresh model length differs")
            sec.L = old_length * float(scale)

            settle_baselines(syn, rows)
            noinput = noinput_local_trace(
                cell, syn, sites_all
            )
            if not np.allclose(
                noinput["t"],
                kinetics["t"],
                rtol=0,
                atol=1e-12,
            ):
                raise RuntimeError(
                    "temporal kinetic/transport grid mismatch"
                )

            local_h, soma_h = impulse_green_matrix(
                cell, syn, sites_all, noinput
            )
            transport = LinearTransport(
                local_h, soma_h
            )

            for timing_name, delays in TIMING_PROGRAMS.items():
                actual = timed_cluster_trace(
                    cell,
                    syn,
                    sites_all,
                    delays,
                    {
                        "t": noinput["t"],
                        "soma": noinput["soma"],
                    },
                )
                ga, gn = timed_conductances(
                    g_ampa_ref,
                    g_nmda_ref,
                    delays,
                )
                reduced = reduced_solve_timed(
                    noinput["local"],
                    transport,
                    ga,
                    gn,
                )
                open_current = synapse_law(
                    noinput["local"],
                    ga,
                    gn,
                )
                open_soma = transport.soma(
                    open_current
                )
                oracle_soma = transport.soma(
                    actual["site_inward_current_nA"]
                )

                current_err = matrix_nrmse(
                    actual["site_inward_current_nA"],
                    reduced["current_nA"],
                )
                reduced_metrics = trace_metrics(
                    actual["soma_depol"],
                    reduced["soma_depol_mV"],
                    actual["t"],
                )
                open_metrics = trace_metrics(
                    actual["soma_depol"],
                    open_soma,
                    actual["t"],
                )
                oracle_metrics = trace_metrics(
                    actual["soma_depol"],
                    oracle_soma,
                    actual["t"],
                )

                converged_flags.append(
                    reduced["converged"]
                )
                spike_flags.append(
                    actual["spike_guard"]
                )

                cases.append(
                    {
                        "geometry_scale": float(scale),
                        "branch_index": int(bi),
                        "section": section_key,
                        "timing_program": timing_name,
                        "delays_ms": list(delays),
                        "reduced_current_nrmse": float(
                            current_err
                        ),
                        "reduced_soma": reduced_metrics,
                        "open_loop_synapse_attacker": open_metrics,
                        "transport_oracle": oracle_metrics,
                        "fixed_point_converged": bool(
                            reduced["converged"]
                        ),
                        "fixed_point_iterations": int(
                            reduced["iterations"]
                        ),
                        "fixed_point_final_relative_update": float(
                            reduced[
                                "final_relative_current_update"
                            ]
                        ),
                        "spike_guard": bool(
                            actual["spike_guard"]
                        ),
                    }
                )

                print(
                    f"g={scale:.2f} [{bi+1}/6] {section_key} "
                    f"{timing_name} "
                    f"J={current_err:.4f} "
                    f"soma={reduced_metrics['nrmse']:.4f} "
                    f"open={open_metrics['nrmse']:.4f} "
                    f"oracle={oracle_metrics['nrmse']:.4f} "
                    f"conv={reduced['converged']}"
                )

    reduced_soma = np.asarray(
        [x["reduced_soma"]["nrmse"] for x in cases],
        dtype=float,
    )
    current_err = np.asarray(
        [x["reduced_current_nrmse"] for x in cases],
        dtype=float,
    )
    oracle_err = np.asarray(
        [x["transport_oracle"]["nrmse"] for x in cases],
        dtype=float,
    )
    open_err = np.asarray(
        [
            x["open_loop_synapse_attacker"]["nrmse"]
            for x in cases
        ],
        dtype=float,
    )

    timing_medians = {}
    for timing_name in TIMING_PROGRAMS:
        vals = [
            x["reduced_soma"]["nrmse"]
            for x in cases
            if x["timing_program"] == timing_name
        ]
        timing_medians[timing_name] = float(
            np.median(vals)
        )

    aggregate = {
        "branches": BRANCHES,
        "geometries": list(GEOMETRIES),
        "timing_programs": {
            k: list(v)
            for k, v in TIMING_PROGRAMS.items()
        },
        "cases": int(len(cases)),
        "multiplicity_per_site": MULTIPLICITY,
        "median_transport_oracle_soma_nrmse": float(
            np.median(oracle_err)
        ),
        "median_reduced_soma_nrmse": float(
            np.median(reduced_soma)
        ),
        "timing_program_median_reduced_soma_nrmse": timing_medians,
        "median_reduced_current_nrmse": float(
            np.median(current_err)
        ),
        "median_open_loop_soma_nrmse": float(
            np.median(open_err)
        ),
        "reduced_to_open_loop_median_error_ratio": float(
            np.median(reduced_soma)
            / (np.median(open_err) + 1e-30)
        ),
        "fraction_reduced_beats_open_loop": float(
            np.mean(reduced_soma < open_err)
        ),
        "fraction_fixed_points_converged": float(
            np.mean(converged_flags)
        ),
        "fraction_spike_guard": float(
            np.mean(spike_flags)
        ),
    }

    transport_ok = (
        aggregate[
            "median_transport_oracle_soma_nrmse"
        ] <= 0.01
    )
    reduced_ok = (
        aggregate["median_reduced_soma_nrmse"] <= 0.02
        and all(
            value <= 0.03
            for value in timing_medians.values()
        )
        and aggregate[
            "median_reduced_current_nrmse"
        ] <= 0.02
        and aggregate[
            "fraction_fixed_points_converged"
        ] == 1.0
        and aggregate["fraction_spike_guard"] == 0.0
    )

    if transport_ok and reduced_ok:
        classification = (
            "TEMPORAL_GREEN_CIRCUIT_GENERALIZES_WITHOUT_REFIT"
        )
        interpretation = (
            "The same Green-matrix plus released-synapse-law reduction "
            "predicts asynchronous forward, reverse and widely staggered "
            "three-site inputs across original and held-out branch metrics "
            "without refitting conductance kinetics or transport."
        )
    elif transport_ok:
        classification = (
            "TRANSPORT_TEMPORALLY_VALID_LOCAL_NONLINEAR_REDUCTION_NOT_PORTABLE"
        )
        interpretation = (
            "The linear transport oracle remains accurate for asynchronous "
            "inputs, but the reduced nonlinear synapse-feedback circuit does "
            "not meet the locked temporal accuracy."
        )
    else:
        classification = (
            "TEMPORAL_GREEN_TRANSPORT_INADEQUATE"
        )
        interpretation = (
            "The geometry-specific Green/transport kernels themselves fail "
            "the temporal oracle control."
        )

    summary = {
        "object": (
            "temporal portability of the local Green-matrix x released "
            "AMPA/NMDA synapse-law reduction"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "geometries": list(GEOMETRIES),
            "timing_programs": {
                k: list(v)
                for k, v in TIMING_PROGRAMS.items()
            },
            "conductance_source": (
                "same single original-geometry HUMAN event used by the "
                "successful 54-case Green-circuit audit; delayed events are "
                "integer-sample shifts only"
            ),
            "fixed_point_damping": DAMPING,
            "fixed_point_tolerance": FP_TOL,
            "fixed_point_max_iterations": FP_MAX_ITER,
            "no_temporal_trace_fit": True,
            "thresholds_locked_before_run": {
                "transport_oracle_median_soma_nrmse_max": 0.01,
                "reduced_median_soma_nrmse_max": 0.02,
                "each_timing_median_soma_nrmse_max": 0.03,
                "reduced_current_median_nrmse_max": 0.02,
                "fixed_point_convergence_fraction_min": 1.0,
                "spike_guard_fraction_max": 0.0,
            },
        },
        "aggregate": aggregate,
        "cases": cases,
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not change delays, analysis window, conductance gain, "
            "damping, geometry scale or temporal alignment after this result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori temporal Green-circuit audit")
    print()
    print(
        "median transport oracle soma NRMSE: "
        f"{aggregate['median_transport_oracle_soma_nrmse']:.4f}"
    )
    print(
        "median reduced soma NRMSE:          "
        f"{aggregate['median_reduced_soma_nrmse']:.4f}"
    )
    print(
        "median reduced current NRMSE:       "
        f"{aggregate['median_reduced_current_nrmse']:.4f}"
    )
    print(
        "median open-loop soma NRMSE:        "
        f"{aggregate['median_open_loop_soma_nrmse']:.4f}"
    )
    print(
        "reduced/open-loop:                  "
        f"{aggregate['reduced_to_open_loop_median_error_ratio']:.4f}"
    )
    print(
        "reduced beats open-loop:            "
        f"{aggregate['fraction_reduced_beats_open_loop']:.3f}"
    )
    for name, value in timing_medians.items():
        print(f"timing {name:12s} median:     {value:.4f}")
    print(
        "fixed-point convergence:            "
        f"{aggregate['fraction_fixed_points_converged']:.3f}"
    )
    print(
        "spike guard:                        "
        f"{aggregate['fraction_spike_guard']:.3f}"
    )
    print(f"classification: {classification}")

    assert len(cases) == 72
    assert np.all(np.isfinite(reduced_soma))
    assert np.all(np.isfinite(current_err))
    assert np.all(np.isfinite(oracle_err))


if __name__ == "__main__":
    main()
