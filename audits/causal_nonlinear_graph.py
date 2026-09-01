from __future__ import annotations

import argparse
import json
import sys
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
from cross_cell_nonlinear_graph import (
    DT_MS,
    EVENT_MS,
    HUMAN_GAMMA,
    N_BRANCHES,
    POST_MS,
    REST_MV,
    TIMING_PROGRAMS,
    first_post_event_indices,
    full_model_trace,
    install_probe_synapses,
    measure_universal_probe_template,
    synapse_law,
    timed_conductances,
)

NEWTON_MAX_ITER = 30
NEWTON_TOL_MV = 1e-10
MAX_BACKTRACK = 8
MIN_ALPHA = 1.0 / 256.0
SITE_CONSISTENCY_TOL_MV = 1e-8


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def nrmse(actual: np.ndarray, pred: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    pred = np.asarray(pred, dtype=float)
    return rms(pred - actual) / (rms(actual) + 1e-30)


def block(v_mV: np.ndarray) -> np.ndarray:
    v = np.asarray(v_mV, dtype=float)
    return 1.0 / (
        1.0 + np.exp(-HUMAN_GAMMA * v) / 3.57
    )


def current_derivative(
    voltage_mV: np.ndarray,
    g_ampa_uS: np.ndarray,
    g_nmda_raw_uS: np.ndarray,
) -> np.ndarray:
    v = np.asarray(voltage_mV, dtype=float)
    ga = np.asarray(g_ampa_uS, dtype=float)
    gn = np.asarray(g_nmda_raw_uS, dtype=float)
    b = block(v)
    bp = HUMAN_GAMMA * b * (1.0 - b)
    # J = [ga + gn*b(v)] * (0-v)
    return gn * bp * (-v) - (ga + gn * b)


def prepare_cell_solver(
    graph: dict,
    cell,
) -> dict:
    G = graph["G_uS"]
    C = np.asarray(graph["C_nF"], dtype=float)
    n = len(C)
    d = C / DT_MS
    A = (
        G
        + sparse.diags(
            d, offsets=0, shape=(n, n), format="csr"
        )
    ).tocsc()
    lu = splu(A)
    soma_node = node_for_section_x(
        graph, cell.soma[0], 0.5
    )
    return {
        "d": d,
        "lu": lu,
        "soma_node": int(soma_node),
        "ncomp": int(n),
    }


def prepare_branch_solver(
    graph: dict,
    cell_solver: dict,
    branch: dict,
) -> dict:
    n = int(cell_solver["ncomp"])
    sec = branch["sec"]
    site_nodes = np.asarray(
        [
            node_for_section_x(
                graph, sec, float(x)
            )
            for x in SITE_X
        ],
        dtype=int,
    )

    B = np.zeros((n, 3), dtype=float)
    for j, node in enumerate(site_nodes):
        B[int(node), int(j)] = 1.0

    X = cell_solver["lu"].solve(B)
    R = X[site_nodes, :]

    return {
        "site_nodes": site_nodes,
        "B": B,
        "X": X,
        "R": R,
    }


def residual_and_jacobian(
    z: np.ndarray,
    passive_sites: np.ndarray,
    R: np.ndarray,
    ga: np.ndarray,
    gn: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    absolute_v = REST_MV + np.asarray(
        z, dtype=float
    )
    current = synapse_law(
        absolute_v,
        ga,
        gn,
    )
    F = z - passive_sites - R @ current
    dJ = current_derivative(
        absolute_v,
        ga,
        gn,
    )
    jac = np.eye(3) - R @ np.diag(dJ)
    return F, jac, current


def solve_site_step(
    passive_sites: np.ndarray,
    R: np.ndarray,
    ga: np.ndarray,
    gn: np.ndarray,
    initial_z: np.ndarray,
) -> dict:
    z = np.asarray(initial_z, dtype=float).copy()
    line_search_failures = 0

    for iteration in range(1, NEWTON_MAX_ITER + 1):
        F, jac, current = residual_and_jacobian(
            z, passive_sites, R, ga, gn
        )
        norm0 = float(
            np.max(np.abs(F))
        )
        if norm0 <= NEWTON_TOL_MV:
            return {
                "z": z,
                "current": current,
                "converged": True,
                "iterations": int(iteration - 1),
                "residual_inf_mV": norm0,
                "line_search_failures": int(
                    line_search_failures
                ),
            }

        try:
            delta = np.linalg.solve(jac, F)
        except np.linalg.LinAlgError:
            return {
                "z": z,
                "current": current,
                "converged": False,
                "iterations": int(iteration),
                "residual_inf_mV": norm0,
                "line_search_failures": int(
                    line_search_failures + 1
                ),
                "failure": "singular_newton_jacobian",
            }

        accepted = False
        alpha = 1.0
        last_candidate = None
        last_norm = None

        for _ in range(MAX_BACKTRACK + 1):
            candidate = z - alpha * delta
            Fc, _, _ = residual_and_jacobian(
                candidate,
                passive_sites,
                R,
                ga,
                gn,
            )
            cand_norm = float(
                np.max(np.abs(Fc))
            )
            last_candidate = candidate
            last_norm = cand_norm
            if cand_norm < norm0:
                z = candidate
                accepted = True
                break
            alpha *= 0.5
            if alpha < MIN_ALPHA - 1e-15:
                break

        if not accepted:
            line_search_failures += 1
            return {
                "z": z,
                "current": current,
                "converged": False,
                "iterations": int(iteration),
                "residual_inf_mV": norm0,
                "candidate_residual_inf_mV": (
                    None
                    if last_norm is None
                    else float(last_norm)
                ),
                "line_search_failures": int(
                    line_search_failures
                ),
                "failure": "backtracking_no_decrease",
            }

    F, _, current = residual_and_jacobian(
        z, passive_sites, R, ga, gn
    )
    return {
        "z": z,
        "current": current,
        "converged": False,
        "iterations": NEWTON_MAX_ITER,
        "residual_inf_mV": float(
            np.max(np.abs(F))
        ),
        "line_search_failures": int(
            line_search_failures
        ),
        "failure": "newton_iteration_limit",
    }


def causal_case(
    cell_solver: dict,
    branch_solver: dict,
    ga: np.ndarray,
    gn: np.ndarray,
    actual_current: np.ndarray,
) -> dict:
    d = cell_solver["d"]
    lu = cell_solver["lu"]
    soma_node = int(cell_solver["soma_node"])

    B = branch_solver["B"]
    X = branch_solver["X"]
    R = branch_solver["R"]
    site_nodes = branch_solver["site_nodes"]

    ncomp = int(cell_solver["ncomp"])
    ntime = int(ga.shape[1])

    state = np.zeros(ncomp, dtype=float)
    oracle_state = np.zeros(ncomp, dtype=float)
    previous_z = np.zeros(3, dtype=float)

    pred_current = np.zeros((3, ntime), dtype=float)
    pred_local_abs = np.zeros((3, ntime), dtype=float)
    pred_soma = np.zeros(ntime, dtype=float)
    oracle_soma = np.zeros(ntime, dtype=float)

    all_converged = True
    max_newton_iterations = 0
    max_residual = 0.0
    total_line_search_failures = 0
    max_site_consistency = 0.0
    first_failure = None

    for ti in range(ntime):
        rhs_pred = d * state
        rhs_oracle = (
            d * oracle_state
            + B @ actual_current[:, ti]
        )
        propagated = lu.solve(
            np.column_stack(
                [rhs_pred, rhs_oracle]
            )
        )
        passive = propagated[:, 0]
        oracle_state = propagated[:, 1]
        passive_sites = passive[site_nodes]

        solved = solve_site_step(
            passive_sites,
            R,
            ga[:, ti],
            gn[:, ti],
            previous_z,
        )

        max_newton_iterations = max(
            max_newton_iterations,
            int(solved["iterations"]),
        )
        max_residual = max(
            max_residual,
            float(solved["residual_inf_mV"]),
        )
        total_line_search_failures += int(
            solved.get("line_search_failures", 0)
        )

        if not solved["converged"]:
            all_converged = False
            if first_failure is None:
                first_failure = {
                    "time_index": int(ti),
                    "time_ms_after_reference": float(
                        (ti + 1) * DT_MS
                    ),
                    "failure": solved.get(
                        "failure", "unknown"
                    ),
                    "residual_inf_mV": float(
                        solved["residual_inf_mV"]
                    ),
                }

        z = np.asarray(solved["z"], dtype=float)
        current = synapse_law(
            REST_MV + z,
            ga[:, ti],
            gn[:, ti],
        )

        state = passive + X @ current
        consistency = float(
            np.max(
                np.abs(
                    state[site_nodes] - z
                )
            )
        )
        max_site_consistency = max(
            max_site_consistency,
            consistency,
        )

        pred_current[:, ti] = current
        pred_local_abs[:, ti] = (
            REST_MV + state[site_nodes]
        )
        pred_soma[ti] = state[soma_node]
        oracle_soma[ti] = oracle_state[soma_node]

        previous_z = state[site_nodes].copy()

    return {
        "current_nA": pred_current,
        "local_voltage_mV": pred_local_abs,
        "soma_depol_mV": pred_soma,
        "oracle_soma_depol_mV": oracle_soma,
        "all_steps_converged": bool(
            all_converged
        ),
        "max_newton_iterations": int(
            max_newton_iterations
        ),
        "max_residual_inf_mV": float(
            max_residual
        ),
        "total_line_search_failures": int(
            total_line_search_failures
        ),
        "max_site_consistency_mV": float(
            max_site_consistency
        ),
        "first_failure": first_failure,
    }


def cell_summary(cases: list[dict]) -> dict:
    out = {}
    for order in sorted(
        {int(x["cell_order"]) for x in cases}
    ):
        rows = [
            x for x in cases
            if int(x["cell_order"]) == order
        ]
        out[str(order)] = {
            "species": rows[0]["species"],
            "morphology_identifier": rows[0][
                "morphology_identifier"
            ],
            "cases": len(rows),
            "median_soma_nrmse": float(
                np.median(
                    [x["causal_soma_nrmse"] for x in rows]
                )
            ),
            "median_current_nrmse": float(
                np.median(
                    [
                        x["causal_current_nrmse"]
                        for x in rows
                    ]
                )
            ),
            "median_oracle_soma_nrmse": float(
                np.median(
                    [
                        x["causal_graph_oracle_nrmse"]
                        for x in rows
                    ]
                )
            ),
            "all_cases_converged": bool(
                all(
                    x["all_steps_converged"]
                    for x in rows
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
            "causal_nonlinear_graph_diagnosis.json"
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
        cell_solver = prepare_cell_solver(
            graph, cell
        )
        branch_solvers = {
            int(branch["branch_index"]):
            prepare_branch_solver(
                graph, cell_solver, branch
            )
            for branch in branches
        }
        probes = install_probe_synapses(branches)

        for branch in branches:
            bi = int(branch["branch_index"])
            solver = branch_solvers[bi]

            for timing_name, delays in TIMING_PROGRAMS.items():
                ga, gn = timed_conductances(
                    g_ampa_ref,
                    g_nmda_ref,
                    delays,
                )
                actual = full_model_trace(
                    cell,
                    probes[bi],
                    delays,
                )
                causal = causal_case(
                    cell_solver,
                    solver,
                    ga,
                    gn,
                    actual["site_inward_current_nA"],
                )

                soma_err = nrmse(
                    actual["soma_depol_mV"],
                    causal["soma_depol_mV"],
                )
                current_err = nrmse(
                    actual["site_inward_current_nA"],
                    causal["current_nA"],
                )
                local_err = nrmse(
                    actual["local_voltage_mV"] - REST_MV,
                    causal["local_voltage_mV"] - REST_MV,
                )
                oracle_err = nrmse(
                    actual["soma_depol_mV"],
                    causal["oracle_soma_depol_mV"],
                )

                cases.append(
                    {
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
                        "causal_soma_nrmse": float(
                            soma_err
                        ),
                        "causal_current_nrmse": float(
                            current_err
                        ),
                        "causal_local_voltage_nrmse": float(
                            local_err
                        ),
                        "causal_graph_oracle_nrmse": float(
                            oracle_err
                        ),
                        "all_steps_converged": bool(
                            causal["all_steps_converged"]
                        ),
                        "max_newton_iterations": int(
                            causal["max_newton_iterations"]
                        ),
                        "max_residual_inf_mV": float(
                            causal["max_residual_inf_mV"]
                        ),
                        "total_line_search_failures": int(
                            causal[
                                "total_line_search_failures"
                            ]
                        ),
                        "max_site_consistency_mV": float(
                            causal[
                                "max_site_consistency_mV"
                            ]
                        ),
                        "first_failure": causal[
                            "first_failure"
                        ],
                        "soma_excursion_guard": bool(
                            actual["excursion_guard"]
                        ),
                    }
                )

                print(
                    f"[{ci+1:02d}/24] "
                    f"{model['species']:5s} "
                    f"{model['morphology_identifier']:>12s} "
                    f"b{bi} {timing_name:12s} "
                    f"soma={soma_err:.4f} "
                    f"J={current_err:.4f} "
                    f"oracle={oracle_err:.4f} "
                    f"conv={causal['all_steps_converged']} "
                    f"it={causal['max_newton_iterations']}"
                )

    soma = np.asarray(
        [x["causal_soma_nrmse"] for x in cases],
        dtype=float,
    )
    current = np.asarray(
        [x["causal_current_nrmse"] for x in cases],
        dtype=float,
    )
    local = np.asarray(
        [x["causal_local_voltage_nrmse"] for x in cases],
        dtype=float,
    )
    oracle = np.asarray(
        [x["causal_graph_oracle_nrmse"] for x in cases],
        dtype=float,
    )
    converged = np.asarray(
        [x["all_steps_converged"] for x in cases],
        dtype=bool,
    )
    guards = np.asarray(
        [x["soma_excursion_guard"] for x in cases],
        dtype=bool,
    )

    per_cell = cell_summary(cases)
    cell_soma = np.asarray(
        [
            x["median_soma_nrmse"]
            for x in per_cell.values()
        ],
        dtype=float,
    )

    timing_medians = {}
    for timing in TIMING_PROGRAMS:
        timing_medians[timing] = float(
            np.median(
                [
                    x["causal_soma_nrmse"]
                    for x in cases
                    if x["timing_program"] == timing
                ]
            )
        )

    species_medians = {}
    for species in ("rat", "human"):
        species_medians[species] = float(
            np.median(
                [
                    x["median_soma_nrmse"]
                    for x in per_cell.values()
                    if x["species"] == species
                ]
            )
        )

    aggregate = {
        "cases": len(cases),
        "median_causal_graph_oracle_soma_nrmse": float(
            np.median(oracle)
        ),
        "median_causal_soma_nrmse": float(
            np.median(soma)
        ),
        "median_causal_current_nrmse": float(
            np.median(current)
        ),
        "median_causal_local_voltage_nrmse": float(
            np.median(local)
        ),
        "timing_median_causal_soma_nrmse": (
            timing_medians
        ),
        "median_cell_causal_soma_nrmse": float(
            np.median(cell_soma)
        ),
        "cells_median_soma_nrmse_le_0p10": int(
            np.sum(cell_soma <= 0.10)
        ),
        "fraction_cases_all_steps_converged": float(
            np.mean(converged)
        ),
        "fraction_soma_excursion_guard": float(
            np.mean(guards)
        ),
        "max_site_consistency_mV": float(
            np.max(
                [
                    x["max_site_consistency_mV"]
                    for x in cases
                ]
            )
        ),
        "max_newton_iterations_any_step": int(
            np.max(
                [
                    x["max_newton_iterations"]
                    for x in cases
                ]
            )
        ),
        "total_line_search_failures": int(
            np.sum(
                [
                    x["total_line_search_failures"]
                    for x in cases
                ]
            )
        ),
        "species_median_cell_soma_nrmse": (
            species_medians
        ),
        "original_global_picard_converged_cases": 207,
        "original_global_picard_total_cases": 288,
    }

    oracle_ok = (
        aggregate[
            "median_causal_graph_oracle_soma_nrmse"
        ] <= 0.02
    )
    causal_ok = (
        aggregate["median_causal_soma_nrmse"] <= 0.03
        and aggregate[
            "median_causal_current_nrmse"
        ] <= 0.03
        and aggregate[
            "median_causal_local_voltage_nrmse"
        ] <= 0.03
        and all(
            value <= 0.05
            for value in timing_medians.values()
        )
        and aggregate[
            "median_cell_causal_soma_nrmse"
        ] <= 0.04
        and aggregate[
            "cells_median_soma_nrmse_le_0p10"
        ] >= 20
        and aggregate[
            "fraction_cases_all_steps_converged"
        ] == 1.0
        and aggregate[
            "fraction_soma_excursion_guard"
        ] == 0.0
    )

    if oracle_ok and causal_ok:
        classification = (
            "CAUSAL_MORPHOLOGY_GRAPH_"
            "NONLINEAR_CLOSURE_VALID"
        )
        failure_mechanism = (
            "GLOBAL_WAVEFORM_PICARD_WAS_THE_"
            "NONPORTABLE_COMPONENT"
        )
    elif oracle_ok:
        classification = (
            "CROSS_CELL_LOCAL_NONLINEAR_"
            "CLOSURE_ITSELF_FAILS"
        )
        failure_mechanism = None
    else:
        classification = (
            "CROSS_CELL_PASSIVE_GRAPH_FAILS_"
            "UNDER_CAUSAL_DRIVE"
        )
        failure_mechanism = None

    strong_label = (
        "SUBPERCENT_CAUSAL_MORPHOLOGY_TO_NONLINEAR_RESPONSE"
        if (
            oracle_ok
            and causal_ok
            and aggregate[
                "median_causal_soma_nrmse"
            ] <= 0.01
        )
        else None
    )

    worst_key = max(
        per_cell,
        key=lambda key: per_cell[key][
            "median_soma_nrmse"
        ],
    )

    summary = {
        "object": (
            "causal state-space diagnosis of cross-cell morphology-graph "
            "plus fixed HUMAN nonlinear probe law"
        ),
        "fci_commit": FCI_COMMIT,
        "protocol": {
            "same_scientific_panel_as_failed_global_picard": True,
            "cells": 24,
            "branches_per_cell": N_BRANCHES,
            "timing_programs": {
                k: list(v)
                for k, v in TIMING_PROGRAMS.items()
            },
            "dt_ms": DT_MS,
            "newton": {
                "max_iterations": NEWTON_MAX_ITER,
                "absolute_inf_residual_tol_mV": NEWTON_TOL_MV,
                "max_backtracking_halvings": MAX_BACKTRACK,
                "minimum_alpha": MIN_ALPHA,
                "site_consistency_tol_mV": SITE_CONSISTENCY_TOL_MV,
                "initial_guess": (
                    "previous solved site depolarization; zero initially"
                ),
            },
            "no_scientific_parameter_changed": True,
            "thresholds_locked_before_run": {
                "oracle_median_nrmse_max": 0.02,
                "soma_median_nrmse_max": 0.03,
                "current_median_nrmse_max": 0.03,
                "local_voltage_median_nrmse_max": 0.03,
                "each_timing_soma_median_nrmse_max": 0.05,
                "median_cell_soma_nrmse_max": 0.04,
                "cells_under_0p10_min": 20,
                "all_steps_converged_fraction_min": 1.0,
                "soma_excursion_guard_fraction_max": 0.0,
                "strong_subpercent_soma_median_max": 0.01,
            },
        },
        "aggregate": aggregate,
        "per_cell": per_cell,
        "worst_cell": per_cell[worst_key],
        "cases": cases,
        "classification": classification,
        "failure_mechanism": failure_mechanism,
        "strong_descriptive_label": strong_label,
        "stopping_line": (
            "Do not tune probe, graph, dt, cell/branch/timing panel, Newton "
            "tolerances or line search after this result."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("Operaattori causal nonlinear graph diagnosis")
    print()
    print(
        "causal graph oracle NRMSE:       "
        f"{aggregate['median_causal_graph_oracle_soma_nrmse']:.4f}"
    )
    print(
        "causal soma NRMSE:               "
        f"{aggregate['median_causal_soma_nrmse']:.4f}"
    )
    print(
        "causal current NRMSE:            "
        f"{aggregate['median_causal_current_nrmse']:.4f}"
    )
    print(
        "causal local voltage NRMSE:      "
        f"{aggregate['median_causal_local_voltage_nrmse']:.4f}"
    )
    for name, value in timing_medians.items():
        print(
            f"timing {name:12s} median:   "
            f"{value:.4f}"
        )
    print(
        "median cell soma NRMSE:          "
        f"{aggregate['median_cell_causal_soma_nrmse']:.4f}"
    )
    print(
        "cells <= 0.10:                   "
        f"{aggregate['cells_median_soma_nrmse_le_0p10']} / 24"
    )
    print(
        "case convergence:                "
        f"{aggregate['fraction_cases_all_steps_converged']:.3f}"
    )
    print(
        "max Newton iterations:           "
        f"{aggregate['max_newton_iterations_any_step']}"
    )
    print(
        "line-search failures:            "
        f"{aggregate['total_line_search_failures']}"
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
        "worst cell:                      "
        f"{per_cell[worst_key]['species']} "
        f"{per_cell[worst_key]['morphology_identifier']} "
        f"{per_cell[worst_key]['median_soma_nrmse']:.4f}"
    )
    print(f"classification: {classification}")
    if failure_mechanism:
        print(f"failure mechanism: {failure_mechanism}")
    if strong_label:
        print(f"strong descriptive label: {strong_label}")

    assert len(cases) == 288
    assert np.all(np.isfinite(soma))
    assert np.all(np.isfinite(current))


if __name__ == "__main__":
    main()
