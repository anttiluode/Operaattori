from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from operaattori.compartment_match import (
    branch_runs,
    greedy_dispersed_match,
    normalized_offdiagonal_coupling,
    run_id_map,
    select_compact_midpoint_sites,
    select_even_sites,
)
from operaattori.nmda_branch import HUMAN, HUMAN_FROZEN_BLOCK, HYBRID_B, solve_equilibrium
from operaattori.real_scaffold import load_morphio_tree
from operaattori.tree_cable import (
    driving_point_impedance_mohm,
    green_impedance_mohm,
    solve_tree_frequency,
)


SOURCE_COMMIT = "75ad8b4d81a7f51bf888b30650c543592340db06"
SOURCE_NAME = "2013_03_06_cell11_1125_H41_06.asc"
SOURCE_REL = (
    "simulating_neurons/neuron_models/human/eyal/"
    "Human_L23_PC_0603_11_937_Eyal_passive_dends_simple_soma/"
    "morphologies/" + SOURCE_NAME
)
SOURCE_URL = (
    "https://raw.githubusercontent.com/ido4848/FCI/"
    + SOURCE_COMMIT + "/" + SOURCE_REL
)


def download_source(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 100_000:
        return dest
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Operaattori-Gate15/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r, dest.open("wb") as f:
        f.write(r.read())
    return dest


def build_electrical_tree(tree):
    n = len(tree.parents)
    parents = tree.parents.copy()
    active = np.zeros(n, dtype=bool)
    active[0] = True
    active[1:] = tree.section_types[1:] != 2

    clamped = np.zeros(n, dtype=bool)
    clamped[0] = True
    for i in range(1, n):
        if active[i] and int(parents[i]) == 0:
            clamped[i] = True

    lengths = np.zeros(n, dtype=float)
    radii = np.zeros(n, dtype=float)
    for i in range(1, n):
        if not active[i] or clamped[i]:
            continue
        p = int(parents[i])
        if p < 0 or not active[p]:
            raise RuntimeError(f"active dendritic node {i} has inactive parent {p}")
        lengths[i] = float(np.linalg.norm(tree.positions[i] - tree.positions[p]))
        radii[i] = max(0.5 * float(tree.radii[i] + tree.radii[p]), 0.15)
        if lengths[i] <= 1e-8:
            active[i] = False

    return parents, lengths, radii, active, clamped


def clamp_current_nA(
    transfer: np.ndarray,
    current_nA: np.ndarray,
) -> float:
    return float(abs(np.dot(np.asarray(transfer, dtype=float), np.asarray(current_nA, dtype=float))))


def isolated_sum_clamp_current(
    Z: np.ndarray,
    transfer: np.ndarray,
    multiplicity: float,
    condition,
) -> float:
    total = 0.0
    for j in range(len(transfer)):
        sol = solve_equilibrium(
            np.asarray([[Z[j, j]]], dtype=float),
            np.asarray([float(multiplicity)]),
            condition,
        )
        if not sol["converged"]:
            raise RuntimeError("isolated-site nonlinear solve failed")
        total += clamp_current_nA(
            np.asarray([transfer[j]]),
            sol["current_nA"],
        )
    return float(total)


def solve_arrangement(
    Z: np.ndarray,
    transfer: np.ndarray,
    multiplicity: float,
    condition,
) -> dict:
    m = np.full(len(transfer), float(multiplicity), dtype=float)
    sol = solve_equilibrium(Z, m, condition)
    current = clamp_current_nA(transfer, sol["current_nA"])
    isolated_sum = isolated_sum_clamp_current(
        Z, transfer, multiplicity, condition
    )
    return {
        "clamp_current_nA": float(current),
        "sum_isolated_clamp_current_nA": float(isolated_sum),
        "interaction_ratio": float(current / (isolated_sum + 1e-30)),
        "max_local_voltage_mV": float(np.max(sol["voltage_mV"])),
        "median_local_voltage_mV": float(np.median(sol["voltage_mV"])),
        "solver_residual_mV": float(sol["residual_mV"]),
        "converged": bool(sol["converged"]),
    }


def plot_locality(out: Path, multiplicities, rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves = np.asarray(
        [[d["locality_index_human_over_frozen"] for d in r["doses"]] for r in rows],
        dtype=float,
    )
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    for curve in curves:
        ax.plot(multiplicities, curve, alpha=0.35)
    ax.plot(
        multiplicities,
        np.median(curves, axis=0),
        marker="o",
        linewidth=2.5,
        label="median branch",
    )
    ax.axhline(1.0, linestyle="--", linewidth=1)
    ax.set_xlabel("synapses per matched site")
    ax.set_ylabel("(human/frozen) clustered / (human/frozen) dispersed")
    ax.set_title("Cell 1125 — NMDA-specific locality index")
    ax.legend()
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--morphology", type=Path)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results" / "gate15")
    ap.add_argument("--branches", type=int, default=10)
    ap.add_argument("--sites", type=int, default=8)
    ap.add_argument("--min-dispersed-runs", type=int, default=6)
    ap.add_argument(
        "--cluster-span-um",
        type=float,
        default=0.0,
        help="if >0, restrict clustered sites to a midpoint window of this physical span",
    )
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if args.sites < 4:
        raise ValueError("Gate 15 needs at least four clustered sites")
    if args.min_dispersed_runs < 3:
        raise ValueError("dispersed control must span at least three branch runs")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(out / "_source" / SOURCE_NAME)
    tree = load_morphio_tree(source)
    parents, lengths, radii, active, clamped = build_electrical_tree(tree)

    state = solve_tree_frequency(
        parents, lengths, radii, active, clamped, 0.0
    )
    runs = branch_runs(parents, state.edge_active, clamped)
    run_ids = run_id_map(len(parents), runs)
    zdrive_complex = driving_point_impedance_mohm(parents, clamped, state)

    if np.max(np.abs(np.imag(zdrive_complex))) > 1e-8:
        raise FloatingPointError("DC driving impedance unexpectedly complex")
    zdrive = np.real(zdrive_complex)
    soma_transfer = np.real(state.transfer_to_clamp)

    pool_sites = np.asarray(
        [
            i for i in range(1, len(parents))
            if state.edge_active[i]
            and run_ids[i] >= 0
            and zdrive[i] > 0
            and abs(soma_transfer[i]) > 0
        ],
        dtype=int,
    )

    candidates = []
    for rid, run in enumerate(runs):
        physical_nodes = [i for i in run if state.edge_active[i]]
        if len(physical_nodes) < args.sites:
            continue
        branch_length = float(np.sum(lengths[physical_nodes]))
        if branch_length < 35.0:
            continue

        if args.cluster_span_um > 0:
            clustered, actual_span = select_compact_midpoint_sites(
                physical_nodes,
                lengths,
                args.sites,
                args.cluster_span_um,
            )
            if len(clustered) != args.sites:
                continue
        else:
            clustered = select_even_sites(physical_nodes, args.sites)
            actual_span = float(branch_length)

        candidates.append(
            (branch_length, rid, physical_nodes, clustered, actual_span)
        )
    candidates.sort(reverse=True, key=lambda x: x[0])

    # Target branches are chosen by geometry only, before seeing nonlinear
    # outcomes. Poor passive matches are reported as a scientific limitation.
    candidates = candidates[: min(args.branches, len(candidates))]
    if len(candidates) < 6:
        raise RuntimeError("too few branches with enough distinct reconstruction nodes")

    multiplicities = np.asarray([1.0, 2.0, 4.0, 6.0])
    conditions = [HUMAN_FROZEN_BLOCK, HYBRID_B, HUMAN]
    rows = []

    for branch_index, (branch_length, rid, run, clustered, actual_cluster_span) in enumerate(candidates):
        dispersed, match_diag = greedy_dispersed_match(
            clustered,
            rid,
            pool_sites,
            run_ids,
            zdrive,
            soma_transfer,
            min_distinct_runs=min(args.min_dispersed_runs, args.sites),
        )

        Zc_complex = green_impedance_mohm(
            parents, clamped, state, clustered
        )
        Zd_complex = green_impedance_mohm(
            parents, clamped, state, dispersed
        )
        complex_ratio = max(
            float(np.max(np.abs(np.imag(Zc_complex))) / (np.max(np.abs(np.real(Zc_complex))) + 1e-30)),
            float(np.max(np.abs(np.imag(Zd_complex))) / (np.max(np.abs(np.real(Zd_complex))) + 1e-30)),
        )
        Zc = np.real(Zc_complex)
        Zd = np.real(Zd_complex)
        Tc = soma_transfer[clustered]
        Td = soma_transfer[dispersed]

        coupling_cluster = normalized_offdiagonal_coupling(Zc)
        coupling_disperse = normalized_offdiagonal_coupling(Zd)

        dose_rows = []
        max_residual = 0.0
        all_converged = True

        for mult in multiplicities:
            solved = {"clustered": {}, "dispersed": {}}
            for cond in conditions:
                solved["clustered"][cond.name] = solve_arrangement(
                    Zc, Tc, float(mult), cond
                )
                solved["dispersed"][cond.name] = solve_arrangement(
                    Zd, Td, float(mult), cond
                )

            for arrangement in ("clustered", "dispersed"):
                for cond in conditions:
                    item = solved[arrangement][cond.name]
                    all_converged &= bool(item["converged"])
                    max_residual = max(max_residual, item["solver_residual_mV"])

            hc = solved["clustered"]["human"]["clamp_current_nA"]
            fc = solved["clustered"]["human_frozen_block"]["clamp_current_nA"]
            hd = solved["dispersed"]["human"]["clamp_current_nA"]
            fd = solved["dispersed"]["human_frozen_block"]["clamp_current_nA"]

            hbc = solved["clustered"]["hybrid_b"]["clamp_current_nA"]
            hbd = solved["dispersed"]["hybrid_b"]["clamp_current_nA"]

            nmda_ratio_cluster = hc / (fc + 1e-30)
            nmda_ratio_disperse = hd / (fd + 1e-30)
            locality = nmda_ratio_cluster / (nmda_ratio_disperse + 1e-30)

            gamma_ratio_cluster = hc / (hbc + 1e-30)
            gamma_ratio_disperse = hd / (hbd + 1e-30)
            gamma_locality = gamma_ratio_cluster / (gamma_ratio_disperse + 1e-30)

            dose_rows.append(
                {
                    "synapses_per_site": float(mult),
                    "total_synapses": int(round(mult * len(clustered))),
                    "human_over_frozen_clustered": float(nmda_ratio_cluster),
                    "human_over_frozen_dispersed": float(nmda_ratio_disperse),
                    "locality_index_human_over_frozen": float(locality),
                    "human_over_hybridB_clustered": float(gamma_ratio_cluster),
                    "human_over_hybridB_dispersed": float(gamma_ratio_disperse),
                    "locality_index_human_over_hybridB": float(gamma_locality),
                    "clustered": solved["clustered"],
                    "dispersed": solved["dispersed"],
                }
            )

        row = {
            "branch_index": int(branch_index),
            "run_id": int(rid),
            "branch_length_um": float(branch_length),
            "cluster_selection_mode": (
                "compact_midpoint_window" if args.cluster_span_um > 0 else "whole_unbranched_run"
            ),
            "requested_cluster_span_um": float(args.cluster_span_um),
            "actual_cluster_span_um": float(actual_cluster_span),
            "clustered_sites": clustered.tolist(),
            "dispersed_sites": dispersed.tolist(),
            "distinct_dispersed_runs": int(match_diag["distinct_match_runs"]),
            "passive_match": match_diag,
            "median_clustered_driving_resistance_Mohm": float(np.median(zdrive[clustered])),
            "median_dispersed_driving_resistance_Mohm": float(np.median(zdrive[dispersed])),
            "median_clustered_abs_soma_transfer": float(np.median(np.abs(Tc))),
            "median_dispersed_abs_soma_transfer": float(np.median(np.abs(Td))),
            "normalized_offdiagonal_coupling_clustered": float(coupling_cluster),
            "normalized_offdiagonal_coupling_dispersed": float(coupling_disperse),
            "coupling_ratio_clustered_over_dispersed": float(
                coupling_cluster / (coupling_disperse + 1e-30)
            ),
            "max_green_imaginary_ratio": float(complex_ratio),
            "all_solvers_converged": bool(all_converged),
            "max_solver_residual_mV": float(max_residual),
            "doses": dose_rows,
        }
        rows.append(row)

        print(
            f"[{branch_index+1:02d}/{len(candidates):02d}] run {rid:3d} "
            f"L={branch_length:6.1f}um span={actual_cluster_span:5.1f}um "
            f"matchZ={match_diag['median_z_ratio_factor']:.3f}x "
            f"matchT={match_diag['median_transfer_ratio_factor']:.3f}x "
            f"runs={match_diag['distinct_match_runs']:2d} "
            f"coupling={row['coupling_ratio_clustered_over_dispersed']:.2f}x "
            f"L48={dose_rows[-1]['locality_index_human_over_frozen']:.4f}"
        )

    high_locality = np.asarray(
        [r["doses"][-1]["locality_index_human_over_frozen"] for r in rows],
        dtype=float,
    )

    def depol(item: dict, field: str = "median_local_voltage_mV") -> float:
        return max(float(item[field]) + 70.0, 1e-12)

    voltage_locality = []
    max_voltage_locality = []
    clustered_human_median_v = []
    clustered_frozen_median_v = []
    dispersed_human_median_v = []
    dispersed_frozen_median_v = []
    for r in rows:
        d = r["doses"][-1]
        ch = d["clustered"]["human"]
        cf = d["clustered"]["human_frozen_block"]
        dh = d["dispersed"]["human"]
        df = d["dispersed"]["human_frozen_block"]

        rc = depol(ch) / depol(cf)
        rd = depol(dh) / depol(df)
        voltage_locality.append(rc / (rd + 1e-30))

        rcm = depol(ch, "max_local_voltage_mV") / depol(cf, "max_local_voltage_mV")
        rdm = depol(dh, "max_local_voltage_mV") / depol(df, "max_local_voltage_mV")
        max_voltage_locality.append(rcm / (rdm + 1e-30))

        clustered_human_median_v.append(ch["median_local_voltage_mV"])
        clustered_frozen_median_v.append(cf["median_local_voltage_mV"])
        dispersed_human_median_v.append(dh["median_local_voltage_mV"])
        dispersed_frozen_median_v.append(df["median_local_voltage_mV"])

    voltage_locality = np.asarray(voltage_locality, dtype=float)
    max_voltage_locality = np.asarray(max_voltage_locality, dtype=float)

    high_gamma_locality = np.asarray(
        [r["doses"][-1]["locality_index_human_over_hybridB"] for r in rows],
        dtype=float,
    )
    coupling_ratio = np.asarray(
        [r["coupling_ratio_clustered_over_dispersed"] for r in rows],
        dtype=float,
    )
    z_match = np.asarray(
        [r["passive_match"]["median_z_ratio_factor"] for r in rows],
        dtype=float,
    )
    t_match = np.asarray(
        [r["passive_match"]["median_transfer_ratio_factor"] for r in rows],
        dtype=float,
    )
    distinct_runs = np.asarray(
        [r["distinct_dispersed_runs"] for r in rows],
        dtype=float,
    )

    dose_median_locality = []
    for di, mult in enumerate(multiplicities):
        vals = np.asarray(
            [r["doses"][di]["locality_index_human_over_frozen"] for r in rows],
            dtype=float,
        )
        dose_median_locality.append(
            {
                "synapses_per_site": float(mult),
                "total_synapses_for_8_sites": int(round(mult * args.sites)),
                "median_locality_index": float(np.median(vals)),
                "fraction_over_1p05": float(np.mean(vals > 1.05)),
            }
        )

    aggregate = {
        "branches": int(len(rows)),
        "sites_per_arrangement": int(args.sites),
        "multiplicities": multiplicities.tolist(),
        "median_passive_z_match_factor": float(np.median(z_match)),
        "median_passive_soma_transfer_match_factor": float(np.median(t_match)),
        "minimum_distinct_dispersed_runs": int(np.min(distinct_runs)),
        "median_coupling_ratio_clustered_over_dispersed": float(np.median(coupling_ratio)),
        "median_high_dose_locality_index_human_over_frozen": float(np.median(high_locality)),
        "median_high_dose_local_voltage_locality_human_over_frozen": float(np.median(voltage_locality)),
        "median_high_dose_max_voltage_locality_human_over_frozen": float(np.median(max_voltage_locality)),
        "median_high_dose_clustered_human_local_voltage_mV": float(np.median(clustered_human_median_v)),
        "median_high_dose_clustered_frozen_local_voltage_mV": float(np.median(clustered_frozen_median_v)),
        "median_high_dose_dispersed_human_local_voltage_mV": float(np.median(dispersed_human_median_v)),
        "median_high_dose_dispersed_frozen_local_voltage_mV": float(np.median(dispersed_frozen_median_v)),
        "fraction_branches_local_voltage_locality_over_1p05": float(np.mean(voltage_locality > 1.05)),
        "fraction_branches_high_dose_locality_over_1p05": float(np.mean(high_locality > 1.05)),
        "fraction_branches_high_dose_locality_over_1p10": float(np.mean(high_locality > 1.10)),
        "median_high_dose_locality_index_human_over_hybridB": float(np.median(high_gamma_locality)),
        "dose_median_locality": dose_median_locality,
        "max_green_imaginary_ratio": float(max(r["max_green_imaginary_ratio"] for r in rows)),
        "max_solver_residual_mV": float(max(r["max_solver_residual_mV"] for r in rows)),
    }

    match_bad = (
        aggregate["median_passive_z_match_factor"] > 1.50
        or aggregate["median_passive_soma_transfer_match_factor"] > 1.50
        or aggregate["minimum_distinct_dispersed_runs"] < min(args.min_dispersed_runs, args.sites)
    )

    if match_bad:
        classification = "PASSIVE_MATCH_INADEQUATE"
        interpretation = (
            "The clustered/dispersed contrast cannot yet be interpreted because the "
            "passive matching or branch dispersion is too loose."
        )
    elif (
        aggregate["median_high_dose_locality_index_human_over_frozen"] >= 1.05
        and aggregate["fraction_branches_high_dose_locality_over_1p05"] >= 0.60
    ):
        classification = "NMDA_LOCALITY_ADVANTAGE_PRESENT"
        interpretation = (
            "After matching individual driving resistance and soma transfer, co-locating "
            "the inputs on one branch amplifies the voltage-dependent NMDA contribution "
            "relative to dispersing the same number of matched inputs across branches."
        )
    else:
        classification = "NO_ROBUST_NMDA_LOCALITY_ADVANTAGE"
        interpretation = (
            "The reduced scaffold does not show a robust extra NMDA benefit for clustered "
            "same-branch inputs after passive site matching."
        )

    summary = {
        "gate": 15,
        "object": "matched clustered-vs-dispersed NMDA compartment audit",
        "source": {
            "paper": "Aizenbud et al. 2026 PNAS",
            "morphology_identifier": "1125",
            "commit": SOURCE_COMMIT,
            "path": SOURCE_REL,
        },
        "assay": {
            "quasi_static_peak_conductance": True,
            "clustered_sites_share_one_maximal_unbranched_run": True,
            "dispersed_sites_exclude_target_run": True,
            "cluster_selection": (
                "compact midpoint physical window"
                if args.cluster_span_um > 0
                else "sites distributed across whole maximal unbranched run"
            ),
            "requested_cluster_span_um": float(args.cluster_span_um),
            "matching_coordinates": [
                "log driving-point impedance",
                "log absolute soma current transfer",
            ],
            "same_synapse_count_and_same_human_conductances": True,
            "frozen_block_preserves_human_AMPA_NMDA_peak_conductances": True,
            "not_a_reproduction_of_FCI": True,
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "branches": rows,
        "stopping_line": (
            "Gate 15 is a mechanistic compartment test, not an FCI estimate. A locality "
            "advantage earns a temporal synaptic assay; failure or poor matching blocks "
            "the jump to developmental growth."
        ),
    }

    (out / "gate15.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    if not args.no_plot:
        plot_locality(
            out / "nmda_locality_index.png",
            multiplicities,
            rows,
        )

    print()
    print("Operaattori Gate 15 — does a branch become a nonlinear compartment?")
    print()
    print(f"branches:                              {aggregate['branches']}")
    print(f"sites / arrangement:                   {aggregate['sites_per_arrangement']}")
    print(f"requested compact span:                {args.cluster_span_um:.1f} um")
    print(f"median passive Z match factor:         {aggregate['median_passive_z_match_factor']:.3f}x")
    print(f"median passive soma-T match factor:    {aggregate['median_passive_soma_transfer_match_factor']:.3f}x")
    print(f"minimum dispersed branch runs:         {aggregate['minimum_distinct_dispersed_runs']}")
    print(f"median clustered/dispersed coupling:   {aggregate['median_coupling_ratio_clustered_over_dispersed']:.3f}x")
    print(f"median high-dose locality H/F:         {aggregate['median_high_dose_locality_index_human_over_frozen']:.4f}")
    print(f"median local-V locality H/F:           {aggregate['median_high_dose_local_voltage_locality_human_over_frozen']:.4f}")
    print(f"median max-V locality H/F:             {aggregate['median_high_dose_max_voltage_locality_human_over_frozen']:.4f}")
    print(f"median clustered V human/frozen:       {aggregate['median_high_dose_clustered_human_local_voltage_mV']:.2f} / {aggregate['median_high_dose_clustered_frozen_local_voltage_mV']:.2f} mV")
    print(f"median dispersed V human/frozen:       {aggregate['median_high_dose_dispersed_human_local_voltage_mV']:.2f} / {aggregate['median_high_dose_dispersed_frozen_local_voltage_mV']:.2f} mV")
    print(f"branches local-V locality >1.05:       {aggregate['fraction_branches_local_voltage_locality_over_1p05']:.3f}")
    print(f"fraction high-dose locality >1.05:     {aggregate['fraction_branches_high_dose_locality_over_1p05']:.3f}")
    print(f"fraction high-dose locality >1.10:     {aggregate['fraction_branches_high_dose_locality_over_1p10']:.3f}")
    print(f"median high-dose locality H/hybrid-B:  {aggregate['median_high_dose_locality_index_human_over_hybridB']:.4f}")
    print(f"max solver residual:                   {aggregate['max_solver_residual_mV']:.3e} mV")
    print()
    print(f"classification: {classification}")
    print(interpretation)

    # Only numerical/protocol integrity assertions. Scientific outcomes may be negative.
    assert len(rows) >= 6
    assert aggregate["max_green_imaginary_ratio"] < 1e-10
    assert aggregate["max_solver_residual_mV"] < 1e-5
    assert np.all(np.isfinite(high_locality))
    assert classification in {
        "PASSIVE_MATCH_INADEQUATE",
        "NMDA_LOCALITY_ADVANTAGE_PRESENT",
        "NO_ROBUST_NMDA_LOCALITY_ADVANTAGE",
    }


if __name__ == "__main__":
    main()
