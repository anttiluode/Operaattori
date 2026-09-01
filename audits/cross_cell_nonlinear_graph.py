from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[1]
AUDITS = ROOT / "audits"
EXPERIMENTS = ROOT / "experiments"
for p in (AUDITS, EXPERIMENTS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from cross_cell_operator import (
    FCI_COMMIT,
    SITE_X,
    choose_branches,
    git_head,
    instantiate_matched_passive,
    load_panel,
    section_name,
    setup_neuron,
)
from direct_cable_graph import (
    build_compartment_graph,
    node_for_section_x,
)
from green_circuit import LinearTransport

DT_MS = 0.025
EVENT_MS = 20.0
POST_MS = 90.0
TSTOP_MS = EVENT_MS + POST_MS + DT_MS
N_BRANCHES = 3

BASE_WEIGHT_US = 0.00088
MULTIPLICITY = 8.0
PROBE_WEIGHT_US = BASE_WEIGHT_US * MULTIPLICITY
HUMAN_RATIO = 0.00131 / 0.00088
HUMAN_GAMMA = 0.078

TAU_R_AMPA = 0.20
TAU_D_AMPA = 1.70
TAU_R_NMDA = 0.29
TAU_D_NMDA = 43.0
REVERSAL_MV = 0.0
REST_MV = -70.0

DAMPING = 0.5
FP_TOL = 1e-8
FP_MAX_ITER = 200
IMPULSE_NA = 0.001

TIMING_PROGRAMS = {
    "synchronous": (0.0, 0.0, 0.0),
    "forward_5": (0.0, 5.0, 10.0),
    "reverse_5": (10.0, 5.0, 0.0),
    "spread_15": (0.0, 15.0, 30.0),
}


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def nrmse(actual: np.ndarray, pred: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    pred = np.asarray(pred, dtype=float)
    return rms(pred - actual) / (rms(actual) + 1e-30)


def block(voltage_mV: np.ndarray) -> np.ndarray:
    v = np.asarray(voltage_mV, dtype=float)
    return 1.0 / (
        1.0 + np.exp(-HUMAN_GAMMA * v) / 3.57
    )


def synapse_law(
    voltage_mV: np.ndarray,
    g_ampa_uS: np.ndarray,
    g_nmda_raw_uS: np.ndarray,
) -> np.ndarray:
    v = np.asarray(voltage_mV, dtype=float)
    ga = np.asarray(g_ampa_uS, dtype=float)
    gn = np.asarray(g_nmda_raw_uS, dtype=float)
    return (
        ga + gn * block(v)
    ) * (REVERSAL_MV - v)


def configure_probe_synapse(syn) -> None:
    syn.e = REVERSAL_MV
    syn.tau_r_AMPA = TAU_R_AMPA
    syn.tau_d_AMPA = TAU_D_AMPA
    syn.tau_r_NMDA = TAU_R_NMDA
    syn.tau_d_NMDA = TAU_D_NMDA
    syn.gamma = HUMAN_GAMMA
    syn.NMDA_ratio = HUMAN_RATIO


def first_post_event_indices(t: np.ndarray) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    expected = int(round(POST_MS / DT_MS))
    idx = np.flatnonzero(
        (t > EVENT_MS + 1e-9)
        & (t <= EVENT_MS + POST_MS + 1e-8)
    )
    if len(idx) != expected:
        raise RuntimeError(
            f"expected {expected} post-event samples, got {len(idx)}; "
            f"first={t[idx[0]] if len(idx) else None}, "
            f"last={t[idx[-1]] if len(idx) else None}"
        )
    rel = t[idx] - EVENT_MS
    expected_rel = (
        np.arange(1, expected + 1, dtype=float) * DT_MS
    )
    if not np.allclose(
        rel, expected_rel, rtol=0, atol=1e-7
    ):
        raise RuntimeError("post-event samples are not on locked fixed grid")
    return idx


def measure_universal_probe_template() -> dict:
    from neuron import h

    sec = h.Section(name="operaattori_probe_template")
    sec.L = 20.0
    sec.diam = 20.0
    sec.nseg = 1
    sec.Ra = 150.0
    sec.cm = 1.0
    sec.insert("pas")
    sec.g_pas = 1.0 / 20000.0
    sec.e_pas = REST_MV

    syn = h.AMPANMDA_EMS(sec(0.5))
    configure_probe_synapse(syn)
    nc = h.NetCon(None, syn)
    nc.weight[0] = PROBE_WEIGHT_US

    # Hold the dummy membrane at rest. Conductance state itself is independent
    # of voltage; the clamp only makes raw-NMDA recovery numerically trivial.
    clamp = h.SEClamp(sec(0.5))
    clamp.dur1 = TSTOP_MS + 10.0
    clamp.amp1 = REST_MV
    clamp.rs = 1e-6

    tvec = h.Vector().record(h._ref_t)
    vvec = h.Vector().record(sec(0.5)._ref_v)
    gavec = h.Vector().record(syn._ref_g_AMPA)
    gnvec = h.Vector().record(syn._ref_g_NMDA)

    h.dt = DT_MS
    h.finitialize(REST_MV)
    h.fcurrent()
    nc.event(EVENT_MS)
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    idx = first_post_event_indices(t)
    v = np.asarray(vvec, dtype=float)[idx]
    ga = np.asarray(gavec, dtype=float)[idx]
    gn_blocked = np.asarray(gnvec, dtype=float)[idx]
    raw_nmda = gn_blocked / np.maximum(block(v), 1e-12)

    receipt = {
        "g_ampa_uS": ga.copy(),
        "g_nmda_raw_uS": raw_nmda.copy(),
        "rest_voltage_max_abs_error_mV": float(
            np.max(np.abs(v - REST_MV))
        ),
        "peak_g_ampa_uS": float(np.max(ga)),
        "peak_g_nmda_raw_uS": float(np.max(raw_nmda)),
    }

    # Drop references before deleting the temporary section.
    del nc, syn, clamp, tvec, vvec, gavec, gnvec
    try:
        h.delete_section(sec=sec)
    except Exception:
        pass

    return receipt


def shift_template(
    template: np.ndarray,
    delay_ms: float,
) -> np.ndarray:
    x = np.asarray(template, dtype=float)
    steps_float = float(delay_ms) / DT_MS
    steps = int(round(steps_float))
    if abs(steps_float - steps) > 1e-9:
        raise ValueError("delay is not an integer number of dt samples")
    out = np.zeros_like(x)
    if steps < len(x):
        out[steps:] = x[: len(x) - steps]
    return out


def timed_conductances(
    g_ampa: np.ndarray,
    g_nmda_raw: np.ndarray,
    delays: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    ga = np.stack(
        [shift_template(g_ampa, d) for d in delays],
        axis=0,
    )
    gn = np.stack(
        [shift_template(g_nmda_raw, d) for d in delays],
        axis=0,
    )
    return ga, gn


def graph_kernels(
    graph: dict,
    cell,
    branches: list[dict],
) -> dict:
    Gmat = graph["G_uS"]
    C = np.asarray(graph["C_nF"], dtype=float)
    ncomp = len(C)

    sources = []
    targets = {}
    for branch in branches:
        bi = int(branch["branch_index"])
        sec = branch["sec"]
        for si, x in enumerate(SITE_X):
            sources.append(
                {
                    "branch_index": bi,
                    "site_index": int(si),
                    "node": node_for_section_x(
                        graph, sec, float(x)
                    ),
                }
            )
            targets[(bi, int(si))] = node_for_section_x(
                graph, sec, float(x)
            )

    ns = len(sources)
    if ns != N_BRANCHES * 3:
        raise RuntimeError(
            f"expected {N_BRANCHES*3} graph sources, got {ns}"
        )

    soma_node = node_for_section_x(
        graph, cell.soma[0], 0.5
    )

    d = C / DT_MS
    A = (
        Gmat
        + sparse.diags(
            d, offsets=0, shape=(ncomp, ncomp), format="csr"
        )
    ).tocsc()
    lu = splu(A)

    ntime = int(round(POST_MS / DT_MS))
    state = np.zeros((ncomp, ns), dtype=float)
    soma = np.zeros((ns, ntime), dtype=float)
    local = {
        key: np.zeros((ns, ntime), dtype=float)
        for key in targets
    }

    impulse = np.zeros((ncomp, ns), dtype=float)
    for col, src in enumerate(sources):
        impulse[int(src["node"]), col] = IMPULSE_NA

    for step in range(ntime):
        rhs = d[:, None] * state
        if step == 0:
            rhs += impulse
        state = lu.solve(rhs)
        soma[:, step] = state[soma_node, :]
        for key, node in targets.items():
            local[key][:, step] = state[node, :]

    packs = {}
    for branch in branches:
        bi = int(branch["branch_index"])
        source_cols = {
            int(src["site_index"]): col
            for col, src in enumerate(sources)
            if int(src["branch_index"]) == bi
        }
        local_h = np.zeros((3, 3, ntime), dtype=float)
        soma_h = np.zeros((3, ntime), dtype=float)
        for sj in range(3):
            col = source_cols[sj]
            soma_h[sj] = soma[col] / IMPULSE_NA
            for ti in range(3):
                local_h[ti, sj] = (
                    local[(bi, ti)][col] / IMPULSE_NA
                )
        packs[bi] = {
            "local_h": local_h,
            "soma_h": soma_h,
        }

    return {
        "packs": packs,
        "ntime": ntime,
    }


def reduced_solve(
    local_h: np.ndarray,
    soma_h: np.ndarray,
    ga: np.ndarray,
    gn: np.ndarray,
) -> dict:
    transport = LinearTransport(local_h, soma_h)
    baseline = np.full_like(ga, REST_MV, dtype=float)

    current = synapse_law(baseline, ga, gn)
    converged = False
    final_error = float("inf")
    iterations = 0

    for iteration in range(1, FP_MAX_ITER + 1):
        voltage = baseline + transport.local(current)
        target = synapse_law(voltage, ga, gn)
        updated = (
            (1.0 - DAMPING) * current
            + DAMPING * target
        )
        final_error = rms(updated - current) / (
            rms(updated) + 1e-30
        )
        current = updated
        iterations = iteration
        if final_error <= FP_TOL:
            converged = True
            break

    voltage = baseline + transport.local(current)
    soma = transport.soma(current)

    return {
        "transport": transport,
        "current_nA": current,
        "local_voltage_mV": voltage,
        "soma_depol_mV": soma,
        "converged": converged,
        "iterations": int(iterations),
        "final_relative_current_update": float(final_error),
    }


def install_probe_synapses(
    branches: list[dict],
) -> dict:
    from neuron import h

    installed = {}
    for branch in branches:
        bi = int(branch["branch_index"])
        sec = branch["sec"]
        syns = []
        ncs = []
        for x in SITE_X:
            syn = h.AMPANMDA_EMS(sec(float(x)))
            configure_probe_synapse(syn)
            nc = h.NetCon(None, syn)
            nc.weight[0] = PROBE_WEIGHT_US
            syns.append(syn)
            ncs.append(nc)
        installed[bi] = {
            "syns": syns,
            "netcons": ncs,
            "sec": sec,
        }
    return installed


def full_model_trace(
    cell,
    probe: dict,
    delays: tuple[float, float, float],
) -> dict:
    from neuron import h

    syns = probe["syns"]
    ncs = probe["netcons"]
    sec = probe["sec"]

    tvec = h.Vector().record(h._ref_t)
    soma_vec = h.Vector().record(cell.soma[0](0.5)._ref_v)
    local_vecs = [
        h.Vector().record(sec(float(x))._ref_v)
        for x in SITE_X
    ]
    i_ampa_vecs = [
        h.Vector().record(syn._ref_i_AMPA)
        for syn in syns
    ]
    i_nmda_vecs = [
        h.Vector().record(syn._ref_i_NMDA)
        for syn in syns
    ]

    h.dt = DT_MS
    h.finitialize(REST_MV)
    h.fcurrent()
    for nc, delay in zip(ncs, delays):
        nc.event(EVENT_MS + float(delay))
    h.continuerun(TSTOP_MS)

    t = np.asarray(tvec, dtype=float)
    idx = first_post_event_indices(t)
    soma_abs = np.asarray(soma_vec, dtype=float)[idx]
    local_abs = np.stack(
        [
            np.asarray(vec, dtype=float)[idx]
            for vec in local_vecs
        ],
        axis=0,
    )
    inward = -(
        np.stack(
            [
                np.asarray(vec, dtype=float)[idx]
                for vec in i_ampa_vecs
            ],
            axis=0,
        )
        + np.stack(
            [
                np.asarray(vec, dtype=float)[idx]
                for vec in i_nmda_vecs
            ],
            axis=0,
        )
    )

    return {
        "soma_depol_mV": soma_abs - REST_MV,
        "soma_abs_mV": soma_abs,
        "local_voltage_mV": local_abs,
        "site_inward_current_nA": inward,
        "excursion_guard": bool(
            np.max(soma_abs) >= -20.0
        ),
        "max_soma_mV": float(np.max(soma_abs)),
    }


def cell_medians(cases: list[dict]) -> dict:
    out = {}
    for cell_order in sorted(
        {int(x["cell_order"]) for x in cases}
    ):
        vals = [
            x for x in cases
            if int(x["cell_order"]) == cell_order
        ]
        out[str(cell_order)] = {
            "species": vals[0]["species"],
            "morphology_identifier": vals[0][
                "morphology_identifier"
            ],
            "cases": len(vals),
            "median_reduced_soma_nrmse": float(
                np.median(
                    [x["reduced_soma_nrmse"] for x in vals]
                )
            ),
            "median_current_nrmse": float(
                np.median(
                    [x["reduced_current_nrmse"] for x in vals]
                )
            ),
            "median_graph_transport_oracle_nrmse": float(
                np.median(
                    [
                        x["graph_transport_oracle_nrmse"]
                        for x in vals
                    ]
                )
            ),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fci-root", type=Path, required=True)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(
            "results/cross_cell_operator/"
            "cross_cell_nonlinear_graph.json"
        ),
    )
    args = ap.parse_args()

    fci_root = args.fci_root.resolve()
    if git_head(fci_root) != FCI_COMMIT:
        raise RuntimeError("FCI source not pinned")

    setup_neuron(fci_root)
    template = measure_universal_probe_template()
    g_ampa_ref = template["g_ampa_uS"]
    g_nmda_ref = template["g_nmda_raw_uS"]

    panel = load_panel(fci_root)
    cases = []
    cell_receipts = []

    for ci, model in enumerate(panel):
        cell = instantiate_matched_passive(
            fci_root, model
        )
        selected = choose_branches(cell)[:N_BRANCHES]
        branches = [
            {
                "branch_index": int(bi),
                "section": section_name(sec),
                "sec": sec,
            }
            for bi, sec in enumerate(selected)
        ]

        graph = build_compartment_graph(cell)
        kernels = graph_kernels(
            graph, cell, branches
        )
        probes = install_probe_synapses(branches)

        cell_case_start = len(cases)

        for branch in branches:
            bi = int(branch["branch_index"])
            pack = kernels["packs"][bi]

            for timing_name, delays in TIMING_PROGRAMS.items():
                ga, gn = timed_conductances(
                    g_ampa_ref,
                    g_nmda_ref,
                    delays,
                )
                reduced = reduced_solve(
                    pack["local_h"],
                    pack["soma_h"],
                    ga,
                    gn,
                )
                actual = full_model_trace(
                    cell,
                    probes[bi],
                    delays,
                )

                current_err = nrmse(
                    actual["site_inward_current_nA"],
                    reduced["current_nA"],
                )
                soma_err = nrmse(
                    actual["soma_depol_mV"],
                    reduced["soma_depol_mV"],
                )

                transport = reduced["transport"]
                oracle_soma = transport.soma(
                    actual["site_inward_current_nA"]
                )
                oracle_err = nrmse(
                    actual["soma_depol_mV"],
                    oracle_soma,
                )

                baseline = np.full_like(
                    ga, REST_MV, dtype=float
                )
                open_current = synapse_law(
                    baseline, ga, gn
                )
                open_soma = transport.soma(
                    open_current
                )
                open_err = nrmse(
                    actual["soma_depol_mV"],
                    open_soma,
                )

                local_err = nrmse(
                    actual["local_voltage_mV"] - REST_MV,
                    reduced["local_voltage_mV"] - REST_MV,
                )

                row = {
                    "cell_order": int(model["order"]),
                    "species": model["species"],
                    "layer": model["layer"],
                    "morphology_identifier": model[
                        "morphology_identifier"
                    ],
                    "branch_index": bi,
                    "section": branch["section"],
                    "timing_program": timing_name,
                    "delays_ms": list(delays),
                    "reduced_soma_nrmse": float(
                        soma_err
                    ),
                    "reduced_current_nrmse": float(
                        current_err
                    ),
                    "reduced_local_voltage_nrmse": float(
                        local_err
                    ),
                    "graph_transport_oracle_nrmse": float(
                        oracle_err
                    ),
                    "open_loop_soma_nrmse": float(
                        open_err
                    ),
                    "reduced_beats_open_loop": bool(
                        soma_err < open_err
                    ),
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
                    "soma_excursion_guard": bool(
                        actual["excursion_guard"]
                    ),
                    "max_soma_mV": actual[
                        "max_soma_mV"
                    ],
                }
                cases.append(row)

                print(
                    f"[{ci+1:02d}/24] "
                    f"{model['species']:5s} "
                    f"{model['morphology_identifier']:>12s} "
                    f"b{bi} {timing_name:12s} "
                    f"soma={soma_err:.4f} "
                    f"J={current_err:.4f} "
                    f"G/T={oracle_err:.4f} "
                    f"open={open_err:.4f} "
                    f"conv={reduced['converged']}"
                )

        cell_rows = cases[cell_case_start:]
        cell_receipts.append(
            {
                "cell_order": int(model["order"]),
                "species": model["species"],
                "layer": model["layer"],
                "morphology_identifier": model[
                    "morphology_identifier"
                ],
                "compartments": int(
                    len(graph["nodes"])
                ),
                "cases": len(cell_rows),
                "median_reduced_soma_nrmse": float(
                    np.median(
                        [
                            x["reduced_soma_nrmse"]
                            for x in cell_rows
                        ]
                    )
                ),
                "median_graph_transport_oracle_nrmse": float(
                    np.median(
                        [
                            x["graph_transport_oracle_nrmse"]
                            for x in cell_rows
                        ]
                    )
                ),
            }
        )

    soma_err = np.asarray(
        [x["reduced_soma_nrmse"] for x in cases],
        dtype=float,
    )
    current_err = np.asarray(
        [x["reduced_current_nrmse"] for x in cases],
        dtype=float,
    )
    local_err = np.asarray(
        [x["reduced_local_voltage_nrmse"] for x in cases],
        dtype=float,
    )
    oracle_err = np.asarray(
        [x["graph_transport_oracle_nrmse"] for x in cases],
        dtype=float,
    )
    open_err = np.asarray(
        [x["open_loop_soma_nrmse"] for x in cases],
        dtype=float,
    )
    converged = np.asarray(
        [x["fixed_point_converged"] for x in cases],
        dtype=bool,
    )
    guards = np.asarray(
        [x["soma_excursion_guard"] for x in cases],
        dtype=bool,
    )

    per_cell = cell_medians(cases)
    cell_soma = np.asarray(
        [
            row["median_reduced_soma_nrmse"]
            for row in per_cell.values()
        ],
        dtype=float,
    )

    timing_medians = {}
    for name in TIMING_PROGRAMS:
        vals = [
            x["reduced_soma_nrmse"]
            for x in cases
            if x["timing_program"] == name
        ]
        timing_medians[name] = float(
            np.median(vals)
        )

    species_medians = {}
    for species in ("rat", "human"):
        vals = [
            row["median_reduced_soma_nrmse"]
            for row in per_cell.values()
            if row["species"] == species
        ]
        species_medians[species] = float(
            np.median(vals)
        )

    aggregate = {
        "cases": len(cases),
        "cells": 24,
        "branches_per_cell": N_BRANCHES,
        "median_graph_transport_oracle_soma_nrmse": float(
            np.median(oracle_err)
        ),
        "median_reduced_soma_nrmse": float(
            np.median(soma_err)
        ),
        "median_reduced_current_nrmse": float(
            np.median(current_err)
        ),
        "median_reduced_local_voltage_nrmse": float(
            np.median(local_err)
        ),
        "median_open_loop_soma_nrmse": float(
            np.median(open_err)
        ),
        "reduced_to_open_loop_median_error_ratio": float(
            np.median(soma_err)
            / (np.median(open_err) + 1e-30)
        ),
        "fraction_reduced_beats_open_loop": float(
            np.mean(soma_err < open_err)
        ),
        "timing_median_reduced_soma_nrmse": (
            timing_medians
        ),
        "median_cell_reduced_soma_nrmse": float(
            np.median(cell_soma)
        ),
        "cells_median_reduced_soma_nrmse_le_0p10": int(
            np.sum(cell_soma <= 0.10)
        ),
        "fraction_fixed_points_converged": float(
            np.mean(converged)
        ),
        "fraction_soma_excursion_guard": float(
            np.mean(guards)
        ),
        "species_median_cell_soma_nrmse": (
            species_medians
        ),
    }

    transport_ok = (
        aggregate[
            "median_graph_transport_oracle_soma_nrmse"
        ] <= 0.02
    )
    reduced_ok = (
        aggregate["median_reduced_soma_nrmse"] <= 0.03
        and aggregate[
            "median_reduced_current_nrmse"
        ] <= 0.03
        and all(
            value <= 0.05
            for value in timing_medians.values()
        )
        and aggregate[
            "median_cell_reduced_soma_nrmse"
        ] <= 0.04
        and aggregate[
            "cells_median_reduced_soma_nrmse_le_0p10"
        ] >= 20
        and aggregate[
            "fraction_fixed_points_converged"
        ] == 1.0
        and aggregate[
            "fraction_reduced_beats_open_loop"
        ] >= 0.90
    )

    if transport_ok and reduced_ok:
        classification = (
            "MORPHOLOGY_GRAPH_X_NONLINEAR_LAW_"
            "PREDICTS_CROSS_CELL_RESPONSE"
        )
    elif transport_ok:
        classification = (
            "CROSS_CELL_GRAPH_TRANSPORT_VALID_"
            "NONLINEAR_CLOSURE_NOT_PORTABLE"
        )
    else:
        classification = (
            "GRAPH_OPERATOR_NOT_ACCURATE_UNDER_"
            "CROSS_CELL_NONLINEAR_DRIVE"
        )

    strong_label = (
        "SUBPERCENT_CROSS_CELL_MORPHOLOGY_TO_NONLINEAR_RESPONSE"
        if (
            transport_ok
            and reduced_ok
            and aggregate[
                "median_reduced_soma_nrmse"
            ] <= 0.01
        )
        else None
    )

    worst_cell_key = max(
        per_cell,
        key=lambda key: per_cell[key][
            "median_reduced_soma_nrmse"
        ],
    )

    summary = {
        "object": (
            "cross-cell nonlinear response prediction from morphology-graph "
            "generated G/T plus a fixed HUMAN AMPA/NMDA probe law"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "cells": 24,
            "branches_per_cell": N_BRANCHES,
            "sites": list(SITE_X),
            "timing_programs": {
                k: list(v)
                for k, v in TIMING_PROGRAMS.items()
            },
            "dt_ms": DT_MS,
            "event_ms": EVENT_MS,
            "post_ms": POST_MS,
            "probe": {
                "label": "HUMAN_PROBE",
                "AMPA_weight_uS": PROBE_WEIGHT_US,
                "base_weight_uS": BASE_WEIGHT_US,
                "multiplicity": MULTIPLICITY,
                "NMDA_ratio": HUMAN_RATIO,
                "gamma_per_mV": HUMAN_GAMMA,
                "tau_r_AMPA_ms": TAU_R_AMPA,
                "tau_d_AMPA_ms": TAU_D_AMPA,
                "tau_r_NMDA_ms": TAU_R_NMDA,
                "tau_d_NMDA_ms": TAU_D_NMDA,
                "reversal_mV": REVERSAL_MV,
            },
            "fixed_point": {
                "damping": DAMPING,
                "relative_tolerance": FP_TOL,
                "max_iterations": FP_MAX_ITER,
            },
            "no_target_electrical_operator_measurement": True,
            "thresholds_locked_before_run": {
                "transport_oracle_median_nrmse_max": 0.02,
                "reduced_soma_median_nrmse_max": 0.03,
                "reduced_current_median_nrmse_max": 0.03,
                "each_timing_median_soma_nrmse_max": 0.05,
                "median_cell_soma_nrmse_max": 0.04,
                "cells_under_0p10_min": 20,
                "fixed_point_convergence_fraction_min": 1.0,
                "reduced_beats_open_loop_fraction_min": 0.90,
                "strong_subpercent_median_soma_nrmse_max": 0.01,
            },
        },
        "universal_conductance_template": {
            "rest_voltage_max_abs_error_mV": template[
                "rest_voltage_max_abs_error_mV"
            ],
            "peak_g_ampa_uS": template[
                "peak_g_ampa_uS"
            ],
            "peak_g_nmda_raw_uS": template[
                "peak_g_nmda_raw_uS"
            ],
        },
        "aggregate": aggregate,
        "per_cell": per_cell,
        "worst_cell": per_cell[worst_cell_key],
        "cell_receipts": cell_receipts,
        "cases": cases,
        "classification": classification,
        "strong_descriptive_label": strong_label,
        "stopping_line": (
            "Do not tune probe weight/gamma/timing, branch selection, passive "
            "constants, graph diameters, damping, alignment, or cell inclusion "
            "after this result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori cross-cell graph x nonlinear-law audit")
    print()
    print(
        "graph transport oracle NRMSE:    "
        f"{aggregate['median_graph_transport_oracle_soma_nrmse']:.4f}"
    )
    print(
        "reduced soma NRMSE:              "
        f"{aggregate['median_reduced_soma_nrmse']:.4f}"
    )
    print(
        "reduced current NRMSE:           "
        f"{aggregate['median_reduced_current_nrmse']:.4f}"
    )
    print(
        "reduced local-voltage NRMSE:     "
        f"{aggregate['median_reduced_local_voltage_nrmse']:.4f}"
    )
    print(
        "open-loop soma NRMSE:            "
        f"{aggregate['median_open_loop_soma_nrmse']:.4f}"
    )
    print(
        "reduced / open-loop:             "
        f"{aggregate['reduced_to_open_loop_median_error_ratio']:.4f}"
    )
    print(
        "reduced beats open-loop:         "
        f"{aggregate['fraction_reduced_beats_open_loop']:.3f}"
    )
    for name, value in timing_medians.items():
        print(
            f"timing {name:12s} median:   "
            f"{value:.4f}"
        )
    print(
        "median cell soma NRMSE:          "
        f"{aggregate['median_cell_reduced_soma_nrmse']:.4f}"
    )
    print(
        "cells <= 0.10:                   "
        f"{aggregate['cells_median_reduced_soma_nrmse_le_0p10']} / 24"
    )
    print(
        "rat median cell soma NRMSE:      "
        f"{species_medians['rat']:.4f}"
    )
    print(
        "human median cell soma NRMSE:    "
        f"{species_medians['human']:.4f}"
    )
    print(
        "fixed-point convergence:         "
        f"{aggregate['fraction_fixed_points_converged']:.3f}"
    )
    print(
        "soma excursion guard:            "
        f"{aggregate['fraction_soma_excursion_guard']:.3f}"
    )
    print(
        "worst cell:                      "
        f"{per_cell[worst_cell_key]['species']} "
        f"{per_cell[worst_cell_key]['morphology_identifier']} "
        f"{per_cell[worst_cell_key]['median_reduced_soma_nrmse']:.4f}"
    )
    print(f"classification: {classification}")
    if strong_label:
        print(f"strong descriptive label: {strong_label}")

    assert len(cases) == 24 * N_BRANCHES * len(TIMING_PROGRAMS)
    assert np.all(np.isfinite(soma_err))
    assert np.all(np.isfinite(current_err))


if __name__ == "__main__":
    main()
