from __future__ import annotations

import argparse
import json
import math
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
    HUMAN_GAMMA,
    dendritic_rows,
    git_head,
)
from gate17_superposition_attack import settle_baselines
from operator_factorization import (
    BRANCHES,
    DT_MS,
    EVENT_MS,
    IMPULSE_NA,
    LENGTH_SCALES,
    MULTIPLICITY,
    POST_MS,
    TSTOP_MS,
    check_branch_identity,
    cluster_trace,
    make_cell,
    recover_branches,
    rms,
    trace_metrics,
)


GEOMETRIES = (1.0, 0.80, 1.20)
PATTERNS = {
    "middle_single": (1,),
    "outer_pair": (0, 2),
    "triple": (0, 1, 2),
}
DAMPING = 0.5
FP_TOL = 1e-8
FP_MAX_ITER = 200


def noinput_local_trace(cell, syn_df, sites: np.ndarray) -> dict:
    from neuron import h

    sites = np.asarray(sites, dtype=int)
    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
    local_vecs = [
        h.Vector().record(
            syn_df.iloc[int(i)]["segments"]._ref_v
        )
        for i in sites
    ]

    h.dt = DT_MS
    h.finitialize(-76.0)
    h.fcurrent()
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    soma = np.asarray(soma_vec, dtype=float)
    local = np.stack(
        [np.asarray(v, dtype=float) for v in local_vecs],
        axis=0,
    )
    post = (t >= EVENT_MS) & (t <= EVENT_MS + POST_MS + 1e-9)
    return {
        "t": t[post] - EVENT_MS,
        "soma": soma[post],
        "local": local[:, post],
    }


def measure_raw_conductance(
    cell,
    syn_df,
    site: int,
) -> dict:
    from neuron import h
    from operator_factorization import configure_human

    site = int(site)
    configure_human(syn_df, np.asarray([site], dtype=int))
    syn = syn_df.iloc[site]["exc_synapses"]
    seg = syn_df.iloc[site]["segments"]

    tvec = h.Vector().record(h._ref_t)
    vvec = h.Vector().record(seg._ref_v)
    gavec = h.Vector().record(syn._ref_g_AMPA)
    gnvec = h.Vector().record(syn._ref_g_NMDA)

    h.dt = DT_MS
    h.finitialize(-76.0)
    h.fcurrent()
    syn_df.iloc[site]["exc_netcons"].event(EVENT_MS)
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    v = np.asarray(vvec, dtype=float)
    ga = np.asarray(gavec, dtype=float)
    gn = np.asarray(gnvec, dtype=float)
    post = (t >= EVENT_MS) & (t <= EVENT_MS + POST_MS + 1e-9)

    vp = v[post]
    block = 1.0 / (
        1.0 + np.exp(-HUMAN_GAMMA * vp) / 3.57
    )
    raw_nmda = gn[post] / np.maximum(block, 1e-12)

    return {
        "t": t[post] - EVENT_MS,
        "g_ampa_uS": ga[post],
        "g_nmda_raw_uS": raw_nmda,
        "reference_voltage_mV": vp,
        "tau": {
            "tau_r_AMPA": float(syn.tau_r_AMPA),
            "tau_d_AMPA": float(syn.tau_d_AMPA),
            "tau_r_NMDA": float(syn.tau_r_NMDA),
            "tau_d_NMDA": float(syn.tau_d_NMDA),
        },
    }


def assert_same_kinetics(
    syn_df,
    sites: np.ndarray,
    reference_tau: dict,
) -> None:
    for i in np.asarray(sites, dtype=int):
        syn = syn_df.iloc[int(i)]["exc_synapses"]
        for name, ref in reference_tau.items():
            value = float(getattr(syn, name))
            if abs(value - ref) > 1e-12:
                raise RuntimeError(
                    f"synaptic kinetic mismatch {name}: {value} != {ref}"
                )


def impulse_green_matrix(
    cell,
    syn_df,
    sites: np.ndarray,
    noinput: dict,
) -> tuple[np.ndarray, np.ndarray]:
    from neuron import h

    sites = np.asarray(sites, dtype=int)
    n = len(noinput["t"])
    local_h = np.zeros((len(sites), len(sites), n), dtype=float)
    soma_h = np.zeros((len(sites), n), dtype=float)

    for sj, source in enumerate(sites):
        source_seg = syn_df.iloc[int(source)]["segments"]
        stim = h.IClamp(
            float(source_seg.x),
            sec=source_seg.sec,
        )
        stim.delay = EVENT_MS
        stim.dur = DT_MS
        stim.amp = IMPULSE_NA

        tvec = h.Vector().record(h._ref_t)
        soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
        local_vecs = [
            h.Vector().record(
                syn_df.iloc[int(i)]["segments"]._ref_v
            )
            for i in sites
        ]

        h.dt = DT_MS
        h.finitialize(-76.0)
        h.fcurrent()
        h.continuerun(TSTOP_MS)

        t = np.asarray(tvec, dtype=float)
        soma = np.asarray(soma_vec, dtype=float)
        local = np.stack(
            [np.asarray(v, dtype=float) for v in local_vecs],
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
            raise RuntimeError("Green kernel time grid differs")

        soma_h[sj] = (
            soma[post] - noinput["soma"]
        ) / IMPULSE_NA
        local_h[:, sj, :] = (
            local[:, post] - noinput["local"]
        ) / IMPULSE_NA

        stim.amp = 0.0

    return local_h, soma_h


def next_pow_two(n: int) -> int:
    return 1 << int(math.ceil(math.log2(max(1, n))))


class LinearTransport:
    def __init__(
        self,
        local_h: np.ndarray,
        soma_h: np.ndarray,
    ):
        self.local_h = np.asarray(local_h, dtype=float)
        self.soma_h = np.asarray(soma_h, dtype=float)
        self.n = int(self.soma_h.shape[1])
        self.nfft = next_pow_two(2 * self.n - 1)
        self.local_f = np.fft.rfft(
            self.local_h,
            n=self.nfft,
            axis=-1,
        )
        self.soma_f = np.fft.rfft(
            self.soma_h,
            n=self.nfft,
            axis=-1,
        )

    def local(self, currents: np.ndarray) -> np.ndarray:
        currents = np.asarray(currents, dtype=float)
        jf = np.fft.rfft(
            currents,
            n=self.nfft,
            axis=-1,
        )
        vf = np.einsum(
            "ijw,jw->iw",
            self.local_f,
            jf,
            optimize=True,
        )
        return np.fft.irfft(
            vf,
            n=self.nfft,
            axis=-1,
        )[:, : self.n]

    def soma(self, currents: np.ndarray) -> np.ndarray:
        currents = np.asarray(currents, dtype=float)
        jf = np.fft.rfft(
            currents,
            n=self.nfft,
            axis=-1,
        )
        vf = np.sum(self.soma_f * jf, axis=0)
        return np.fft.irfft(
            vf,
            n=self.nfft,
        )[: self.n]


def synapse_law(
    voltage_mV: np.ndarray,
    g_ampa_uS: np.ndarray,
    g_nmda_raw_uS: np.ndarray,
) -> np.ndarray:
    block = 1.0 / (
        1.0
        + np.exp(
            -HUMAN_GAMMA
            * np.asarray(voltage_mV, dtype=float)
        )
        / 3.57
    )
    inward = (
        np.asarray(g_ampa_uS, dtype=float)
        + np.asarray(g_nmda_raw_uS, dtype=float) * block
    ) * (-np.asarray(voltage_mV, dtype=float))
    return inward


def reduced_solve(
    noinput_local: np.ndarray,
    transport: LinearTransport,
    g_ampa_ref: np.ndarray,
    g_nmda_ref: np.ndarray,
    active_local_indices: np.ndarray,
) -> dict:
    noinput_local = np.asarray(noinput_local, dtype=float)
    nsite, n = noinput_local.shape
    ga = np.zeros((nsite, n), dtype=float)
    gn = np.zeros((nsite, n), dtype=float)
    for i in np.asarray(active_local_indices, dtype=int):
        ga[int(i)] = g_ampa_ref
        gn[int(i)] = g_nmda_ref

    current = synapse_law(
        noinput_local,
        ga,
        gn,
    )
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


def matrix_nrmse(
    actual: np.ndarray,
    pred: np.ndarray,
) -> float:
    actual = np.asarray(actual, dtype=float)
    pred = np.asarray(pred, dtype=float)
    return rms(pred - actual) / (rms(actual) + 1e-30)


def compact_case(
    geometry: float,
    section: str,
    pattern: str,
    current_nrmse: float,
    reduced_metrics: dict,
    oracle_metrics: dict,
    converged: bool,
    iterations: int,
    final_err: float,
    frozen_soma_metrics: dict | None,
    frozen_current_metrics: dict | None,
) -> dict:
    row = {
        "geometry_scale": float(geometry),
        "section": section,
        "pattern": pattern,
        "reduced_current_nrmse": float(current_nrmse),
        "reduced_soma": reduced_metrics,
        "transport_oracle": oracle_metrics,
        "fixed_point_converged": bool(converged),
        "fixed_point_iterations": int(iterations),
        "fixed_point_final_relative_update": float(final_err),
    }
    if frozen_soma_metrics is not None:
        row["frozen_soma_attacker"] = frozen_soma_metrics
    if frozen_current_metrics is not None:
        row["frozen_current_factorization"] = frozen_current_metrics
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument("--cluster-span-um", type=float, default=55.0)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/operator_factorization/green_circuit.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source is not pinned")

    # Recover branch identities and one universal synapse-kinetics trace from
    # the untouched released model.
    ref_cell, ref_syn = make_cell(fci_root)
    ref_rows, branches = recover_branches(
        ref_cell, ref_syn, args.cluster_span_um
    )
    settle_baselines(ref_syn, ref_rows)
    ref_site = int(branches[0]["sites"][1])
    kinetics = measure_raw_conductance(
        ref_cell,
        ref_syn,
        ref_site,
    )
    g_ampa_ref = np.asarray(
        kinetics["g_ampa_uS"],
        dtype=float,
    )
    g_nmda_ref = np.asarray(
        kinetics["g_nmda_raw_uS"],
        dtype=float,
    )

    for branch in branches:
        assert_same_kinetics(
            ref_syn,
            np.asarray(branch["sites"], dtype=int),
            kinetics["tau"],
        )

    original_lookup = {}
    cases = []
    all_converged = []
    all_spikes = []

    for bi, branch in enumerate(branches):
        section_key = branch["canonical_section"]
        original_lookup[section_key] = {}

        for scale in GEOMETRIES:
            cell, syn = make_cell(fci_root)
            rows = dendritic_rows(syn)
            sites_all = np.asarray(branch["sites"], dtype=int)
            check_branch_identity(
                syn, sites_all, section_key
            )
            assert_same_kinetics(
                syn,
                sites_all,
                kinetics["tau"],
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
                cell,
                syn,
                sites_all,
            )
            if not np.allclose(
                noinput["t"],
                kinetics["t"],
                rtol=0,
                atol=1e-12,
            ):
                raise RuntimeError(
                    "kinetic/transport time grid mismatch"
                )

            local_h, soma_h = impulse_green_matrix(
                cell,
                syn,
                sites_all,
                noinput,
            )
            transport = LinearTransport(
                local_h,
                soma_h,
            )

            for pattern_name, local_indices_tuple in PATTERNS.items():
                local_indices = np.asarray(
                    local_indices_tuple,
                    dtype=int,
                )
                active_sites = sites_all[local_indices]
                actual = cluster_trace(
                    cell,
                    syn,
                    active_sites,
                    {
                        "t": noinput["t"],
                        "soma": noinput["soma"],
                    },
                )
                all_spikes.append(actual["spike_guard"])

                reduced = reduced_solve(
                    noinput["local"],
                    transport,
                    g_ampa_ref,
                    g_nmda_ref,
                    local_indices,
                )
                all_converged.append(reduced["converged"])

                actual_current_full = np.zeros(
                    (3, len(noinput["t"])),
                    dtype=float,
                )
                actual_current_full[local_indices] = (
                    actual["site_inward_current_nA"]
                )
                current_err = matrix_nrmse(
                    actual_current_full,
                    reduced["current_nA"],
                )
                reduced_metrics = trace_metrics(
                    actual["soma_depol"],
                    reduced["soma_depol_mV"],
                    actual["t"],
                )
                oracle_soma = transport.soma(
                    actual_current_full
                )
                oracle_metrics = trace_metrics(
                    actual["soma_depol"],
                    oracle_soma,
                    actual["t"],
                )

                frozen_soma_metrics = None
                frozen_current_metrics = None
                if abs(scale - 1.0) < 1e-12:
                    original_lookup[section_key][pattern_name] = {
                        "soma": actual["soma_depol"].copy(),
                        "current_full": actual_current_full.copy(),
                    }
                else:
                    base = original_lookup[
                        section_key
                    ][pattern_name]
                    frozen_soma_metrics = trace_metrics(
                        actual["soma_depol"],
                        base["soma"],
                        actual["t"],
                    )
                    frozen_current_soma = transport.soma(
                        base["current_full"]
                    )
                    frozen_current_metrics = trace_metrics(
                        actual["soma_depol"],
                        frozen_current_soma,
                        actual["t"],
                    )

                cases.append(
                    compact_case(
                        scale,
                        section_key,
                        pattern_name,
                        current_err,
                        reduced_metrics,
                        oracle_metrics,
                        reduced["converged"],
                        reduced["iterations"],
                        reduced[
                            "final_relative_current_update"
                        ],
                        frozen_soma_metrics,
                        frozen_current_metrics,
                    )
                )

                print(
                    f"g={scale:.2f} [{bi+1}/6] "
                    f"{section_key} {pattern_name} "
                    f"J={current_err:.4f} "
                    f"soma={reduced_metrics['nrmse']:.4f} "
                    f"oracle={oracle_metrics['nrmse']:.4f} "
                    f"conv={reduced['converged']} "
                    f"it={reduced['iterations']}"
                )

    soma_err = np.asarray(
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
    heldout = [
        x for x in cases
        if abs(x["geometry_scale"] - 1.0) > 1e-12
    ]
    frozen_current_err = np.asarray(
        [
            x["frozen_current_factorization"]["nrmse"]
            for x in heldout
        ],
        dtype=float,
    )
    reduced_heldout_err = np.asarray(
        [x["reduced_soma"]["nrmse"] for x in heldout],
        dtype=float,
    )

    pattern_medians = {}
    for pattern_name in PATTERNS:
        vals = [
            x["reduced_soma"]["nrmse"]
            for x in cases
            if x["pattern"] == pattern_name
        ]
        pattern_medians[pattern_name] = float(
            np.median(vals)
        )

    aggregate = {
        "branches": BRANCHES,
        "geometries": list(GEOMETRIES),
        "patterns": list(PATTERNS),
        "cases": int(len(cases)),
        "heldout_cases": int(len(heldout)),
        "median_transport_oracle_soma_nrmse": float(
            np.median(oracle_err)
        ),
        "median_reduced_soma_nrmse": float(
            np.median(soma_err)
        ),
        "pattern_median_reduced_soma_nrmse": pattern_medians,
        "median_reduced_current_nrmse": float(
            np.median(current_err)
        ),
        "median_heldout_frozen_current_factorization_nrmse": float(
            np.median(frozen_current_err)
        ),
        "median_heldout_reduced_soma_nrmse": float(
            np.median(reduced_heldout_err)
        ),
        "reduced_to_frozen_current_median_error_ratio": float(
            np.median(reduced_heldout_err)
            / (np.median(frozen_current_err) + 1e-30)
        ),
        "fraction_reduced_beats_frozen_current": float(
            np.mean(
                reduced_heldout_err < frozen_current_err
            )
        ),
        "fraction_fixed_points_converged": float(
            np.mean(all_converged)
        ),
        "fraction_actual_spike_guard": float(
            np.mean(all_spikes)
        ),
    }

    pass_all = (
        aggregate["median_transport_oracle_soma_nrmse"]
        <= 0.02
        and aggregate["median_reduced_soma_nrmse"] <= 0.05
        and all(
            v <= 0.08
            for v in pattern_medians.values()
        )
        and aggregate["median_reduced_current_nrmse"] <= 0.08
        and aggregate[
            "reduced_to_frozen_current_median_error_ratio"
        ] <= 0.80
        and aggregate[
            "fraction_reduced_beats_frozen_current"
        ] >= 24.0 / 36.0
        and aggregate[
            "fraction_fixed_points_converged"
        ] == 1.0
        and aggregate["fraction_actual_spike_guard"] == 0.0
    )

    if pass_all:
        classification = (
            "LOCAL_GREEN_MATRIX_X_SYNAPSE_LAW_REDUCES_RELEASED_NEURON"
        )
        interpretation = (
            "The released HUMAN synapse law, coupled only through measured "
            "geometry-specific local and soma impulse kernels, predicts the "
            "full model across branches, input patterns and held-out branch "
            "length perturbations without branch-specific nonlinear waveform "
            "lookups."
        )
    elif (
        aggregate["median_transport_oracle_soma_nrmse"]
        <= 0.02
    ):
        classification = (
            "TRANSPORT_REDUCTION_VALID_SYNAPSE_FEEDBACK_REDUCTION_INADEQUATE"
        )
        interpretation = (
            "Linear transport remains accurate, but the self-consistent "
            "three-site synapse-law reduction does not reproduce the full "
            "released model at the preregistered accuracy."
        )
    else:
        classification = (
            "GREEN_TRANSPORT_KERNELS_INADEQUATE"
        )
        interpretation = (
            "The measured linear Green/impulse kernels themselves fail the "
            "oracle transport control."
        )

    summary = {
        "object": (
            "three-site nonlinear conductance circuit coupled by measured "
            "local Green matrix and site-to-soma transport kernels"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "geometries": list(GEOMETRIES),
            "patterns": {
                k: list(v) for k, v in PATTERNS.items()
            },
            "multiplicity_per_active_site": MULTIPLICITY,
            "human_gamma_per_mV": HUMAN_GAMMA,
            "conductance_kinetics_source": (
                "one original-geometry HUMAN synapse event; raw NMDA recovered "
                "by dividing g_NMDA by the released magnesium block"
            ),
            "fixed_point_damping": DAMPING,
            "fixed_point_relative_tolerance": FP_TOL,
            "fixed_point_max_iterations": FP_MAX_ITER,
            "no_branch_specific_nonlinear_fit": True,
            "no_heldout_trace_fit": True,
            "thresholds_locked_before_run": {
                "transport_oracle_median_nrmse_max": 0.02,
                "reduced_soma_median_nrmse_max": 0.05,
                "each_pattern_reduced_soma_median_nrmse_max": 0.08,
                "reduced_current_median_nrmse_max": 0.08,
                "reduced_to_frozen_current_median_error_ratio_max": 0.80,
                "reduced_beats_frozen_current_fraction_min": 24.0 / 36.0,
                "fixed_point_convergence_fraction_min": 1.0,
                "spike_guard_fraction_max": 0.0,
            },
        },
        "kinetics": {
            "tau": kinetics["tau"],
            "samples": int(len(kinetics["t"])),
        },
        "aggregate": aggregate,
        "cases": cases,
        "classification": classification,
        "interpretation": interpretation,
        "stopping_line": (
            "Do not fit conductance gains, change damping, align traces, "
            "select branches or alter the pattern/geometry panel after this "
            "result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori local Green matrix x synapse-law audit")
    print()
    print(
        "median transport oracle soma NRMSE:   "
        f"{aggregate['median_transport_oracle_soma_nrmse']:.4f}"
    )
    print(
        "median reduced soma NRMSE:            "
        f"{aggregate['median_reduced_soma_nrmse']:.4f}"
    )
    print(
        "median reduced current NRMSE:         "
        f"{aggregate['median_reduced_current_nrmse']:.4f}"
    )
    print(
        "heldout frozen-current NRMSE:         "
        f"{aggregate['median_heldout_frozen_current_factorization_nrmse']:.4f}"
    )
    print(
        "heldout reduced soma NRMSE:           "
        f"{aggregate['median_heldout_reduced_soma_nrmse']:.4f}"
    )
    print(
        "reduced / frozen-current:             "
        f"{aggregate['reduced_to_frozen_current_median_error_ratio']:.4f}"
    )
    print(
        "reduced beats frozen-current:         "
        f"{aggregate['fraction_reduced_beats_frozen_current']:.3f}"
    )
    for name, value in pattern_medians.items():
        print(f"pattern {name:16s} median:    {value:.4f}")
    print(
        "fixed-point convergence:              "
        f"{aggregate['fraction_fixed_points_converged']:.3f}"
    )
    print(
        "spike guard:                          "
        f"{aggregate['fraction_actual_spike_guard']:.3f}"
    )
    print(f"classification: {classification}")

    assert len(cases) == 54
    assert np.all(np.isfinite(soma_err))
    assert np.all(np.isfinite(current_err))
    assert np.all(np.isfinite(oracle_err))


if __name__ == "__main__":
    main()
