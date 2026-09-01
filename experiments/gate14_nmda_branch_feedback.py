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

from operaattori.nmda_branch import (
    HUMAN, HYBRID_B, HUMAN_FROZEN_BLOCK, HUMAN_LINEAR_CURRENT,
    effective_rank_rows, linear_regression_r2, solve_equilibrium,
)
from operaattori.real_scaffold import load_morphio_tree
from operaattori.tree_cable import (
    active_child_counts, green_impedance_mohm, solve_tree_frequency,
)


SOURCE_COMMIT = "75ad8b4d81a7f51bf888b30650c543592340db06"
SOURCE_NAME = "2013_03_06_cell11_1125_H41_06.asc"
SOURCE_REL = (
    "simulating_neurons/neuron_models/human/eyal/"
    "Human_L23_PC_0603_11_937_Eyal_passive_dends_simple_soma/"
    "morphologies/" + SOURCE_NAME
)
SOURCE_URL = "https://raw.githubusercontent.com/ido4848/FCI/" + SOURCE_COMMIT + "/" + SOURCE_REL


def download_source(dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 100_000:
        return dest
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Operaattori-Gate14/1.0"})
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
            raise RuntimeError(f'active node {i} has inactive parent {p}')
        lengths[i] = float(np.linalg.norm(tree.positions[i] - tree.positions[p]))
        radii[i] = max(0.5 * float(tree.radii[i] + tree.radii[p]), 0.15)
        if lengths[i] <= 1e-8:
            active[i] = False
    return parents, lengths, radii, active, clamped


def branch_runs(parents, state, clamped):
    children: dict[int, list[int]] = {}
    for i in range(1, len(parents)):
        if state.edge_active[i]:
            children.setdefault(int(parents[i]), []).append(i)

    starts = set(int(i) for i in np.flatnonzero(clamped))
    starts.update(
        i for i in children
        if not clamped[i] and len(children.get(i, [])) != 1
    )

    runs = []
    for start in sorted(starts):
        for child in children.get(start, []):
            path = [child]
            cur = child
            while len(children.get(cur, [])) == 1:
                cur = children[cur][0]
                path.append(cur)
            runs.append(path)
    return runs


def synapse_slots(path, lengths):
    slots = []
    for node in path:
        count = max(1, int(round(float(lengths[node]))))
        slots.extend([int(node)] * count)
    return np.asarray(slots, dtype=int)


def select_slots(slots: np.ndarray, n: int) -> np.ndarray:
    n = min(int(n), len(slots))
    if n <= 0:
        return np.empty(0, dtype=int)
    idx = np.rint(np.linspace(0, len(slots) - 1, n)).astype(int)
    return slots[idx]


def multiplicity_on_union(chosen: np.ndarray, union_sites: np.ndarray) -> np.ndarray:
    lookup = {int(s): i for i, s in enumerate(union_sites)}
    m = np.zeros(len(union_sites), dtype=float)
    for s in chosen:
        m[lookup[int(s)]] += 1.0
    return m


def plot_curves(path: Path, doses, condition_curves, branch_label: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    for name, curve in condition_curves.items():
        ax.plot(doses, curve, marker='o', label=name)
    ax.set_xlabel('simultaneous synapses distributed on one branch')
    ax.set_ylabel('|soma-clamp current| (nA)')
    ax.set_title(branch_label)
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--morphology', type=Path)
    ap.add_argument('--out-dir', type=Path, default=ROOT/'results'/'gate14')
    ap.add_argument('--branches', type=int, default=12)
    ap.add_argument('--no-plot', action='store_true')
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(out/'_source'/SOURCE_NAME)
    tree = load_morphio_tree(source)
    parents, lengths, radii, active, clamped = build_electrical_tree(tree)
    state = solve_tree_frequency(parents, lengths, radii, active, clamped, 0.0)

    runs = branch_runs(parents, state, clamped)
    candidates = []
    for path in runs:
        L = float(np.sum(lengths[path]))
        slots = synapse_slots(path, lengths)
        if L >= 40.0 and len(slots) >= 48:
            candidates.append((L, path, slots))
    candidates.sort(key=lambda x: x[0], reverse=True)
    candidates = candidates[:min(args.branches, len(candidates))]
    if len(candidates) < 6:
        raise RuntimeError('too few long branches for Gate 14')

    doses = np.asarray([1, 2, 4, 8, 16, 32, 48], dtype=int)
    conditions = [HUMAN_LINEAR_CURRENT, HUMAN_FROZEN_BLOCK, HYBRID_B, HUMAN]
    rows = []
    curves_by_condition = {c.name: [] for c in conditions}
    most_effect = None

    for bi, (branch_length, path, slots) in enumerate(candidates):
        selections = {int(n): select_slots(slots, int(n)) for n in doses}
        union_sites = np.unique(np.concatenate(list(selections.values())))
        Zc = green_impedance_mohm(parents, clamped, state, union_sites)
        imag_ratio = float(np.max(np.abs(np.imag(Zc))) / (np.max(np.abs(np.real(Zc))) + 1e-30))
        Z = np.real(Zc)
        symmetry_error = float(np.max(np.abs(Z - Z.T)))
        T = np.real(state.transfer_to_clamp[union_sites])

        # One-synapse isolated response per site and condition, used as a
        # superposition/cooperativity ruler.
        single_clamp = {}
        for cond in conditions:
            values = np.zeros(len(union_sites), dtype=float)
            for j in range(len(union_sites)):
                sol = solve_equilibrium(
                    np.asarray([[Z[j,j]]]), np.asarray([1.0]), cond
                )
                if not sol['converged']:
                    raise RuntimeError(f'single-site solve failed on branch {bi}')
                values[j] = abs(T[j] * sol['current_nA'][0])
            single_clamp[cond.name] = values

        condition_curves = {c.name: [] for c in conditions}
        condition_coop = {c.name: [] for c in conditions}
        condition_max_v = {c.name: [] for c in conditions}
        all_converged = True
        max_residual = 0.0

        for n in doses:
            mult = multiplicity_on_union(selections[int(n)], union_sites)
            active_idx = np.flatnonzero(mult > 0)
            Za = Z[np.ix_(active_idx, active_idx)]
            Ta = T[active_idx]
            ma = mult[active_idx]

            for cond in conditions:
                sol = solve_equilibrium(Za, ma, cond)
                all_converged &= bool(sol['converged'])
                max_residual = max(max_residual, float(sol['residual_mV']))
                clamp = float(abs(np.dot(Ta, sol['current_nA'])))
                independent = float(np.sum(ma * single_clamp[cond.name][active_idx]))
                coop = clamp / (independent + 1e-30)
                condition_curves[cond.name].append(clamp)
                condition_coop[cond.name].append(coop)
                condition_max_v[cond.name].append(float(np.max(sol['voltage_mV'])))

        for cond in conditions:
            curves_by_condition[cond.name].append(condition_curves[cond.name])

        human = np.asarray(condition_curves['human'])
        frozen = np.asarray(condition_curves['human_frozen_block'])
        hybrid = np.asarray(condition_curves['hybrid_b'])
        linear = np.asarray(condition_curves['human_linear_current'])
        human_coop = np.asarray(condition_coop['human'])
        frozen_coop = np.asarray(condition_coop['human_frozen_block'])
        linear_coop = np.asarray(condition_coop['human_linear_current'])

        diagZ = np.diag(Z)
        row = {
            'branch_index': int(bi),
            'proximal_node': int(parents[path[0]]),
            'distal_node': int(path[-1]),
            'branch_length_um': float(branch_length),
            'path_nodes': int(len(path)),
            'synapse_slots_approx_1_per_um': int(len(slots)),
            'union_sites': int(len(union_sites)),
            'median_driving_point_resistance_Mohm': float(np.median(diagZ)),
            'max_driving_point_resistance_Mohm': float(np.max(diagZ)),
            'green_symmetry_error_Mohm': symmetry_error,
            'green_imaginary_ratio': imag_ratio,
            'all_solvers_converged': bool(all_converged),
            'max_solver_residual_mV': float(max_residual),
            'human_to_frozen_high_dose_ratio': float(human[-1]/(frozen[-1]+1e-30)),
            'human_to_hybridB_high_dose_ratio': float(human[-1]/(hybrid[-1]+1e-30)),
            'human_high_dose_cooperativity': float(human_coop[-1]),
            'frozen_high_dose_cooperativity': float(frozen_coop[-1]),
            'linear_high_dose_cooperativity': float(linear_coop[-1]),
            'human_peak_local_voltage_mV': float(np.max(condition_max_v['human'])),
            'curves': {name: [float(x) for x in values] for name, values in condition_curves.items()},
            'cooperativity': {name: [float(x) for x in values] for name, values in condition_coop.items()},
        }
        rows.append(row)
        effect = row['human_to_frozen_high_dose_ratio']
        if most_effect is None or effect > most_effect[0]:
            most_effect = (effect, bi, condition_curves)

    for name in curves_by_condition:
        curves_by_condition[name] = np.asarray(curves_by_condition[name], dtype=float)

    human_enh = np.asarray([r['human_to_frozen_high_dose_ratio'] for r in rows])
    human_hybrid = np.asarray([r['human_to_hybridB_high_dose_ratio'] for r in rows])
    human_coop = np.asarray([r['human_high_dose_cooperativity'] for r in rows])
    frozen_coop = np.asarray([r['frozen_high_dose_cooperativity'] for r in rows])
    linear_coop = np.asarray([r['linear_high_dose_cooperativity'] for r in rows])
    Rmed = np.asarray([r['median_driving_point_resistance_Mohm'] for r in rows])
    Rmax = np.asarray([r['max_driving_point_resistance_Mohm'] for r in rows])
    L = np.asarray([r['branch_length_um'] for r in rows])

    ranks = {name: effective_rank_rows(curve) for name, curve in curves_by_condition.items()}
    scalar_r2 = linear_regression_r2(np.column_stack([Rmed, Rmax, L]), human_enh)

    aggregate = {
        'branches': int(len(rows)),
        'doses': doses.tolist(),
        'median_human_to_frozen_high_dose_ratio': float(np.median(human_enh)),
        'max_human_to_frozen_high_dose_ratio': float(np.max(human_enh)),
        'median_human_to_hybridB_high_dose_ratio': float(np.median(human_hybrid)),
        'fraction_branches_human_over_frozen_10pct': float(np.mean(human_enh > 1.10)),
        'median_human_high_dose_cooperativity': float(np.median(human_coop)),
        'median_frozen_high_dose_cooperativity': float(np.median(frozen_coop)),
        'median_human_minus_frozen_cooperativity': float(np.median(human_coop - frozen_coop)),
        'max_abs_linear_cooperativity_minus_one': float(np.max(np.abs(linear_coop - 1.0))),
        'dose_curve_effective_rank': {k: float(v) for k,v in ranks.items()},
        'human_minus_frozen_effective_rank': float(ranks['human'] - ranks['human_frozen_block']),
        'human_enhancement_r2_from_Rmed_Rmax_length': float(scalar_r2),
        'max_green_symmetry_error_Mohm': float(max(r['green_symmetry_error_Mohm'] for r in rows)),
        'max_green_imaginary_ratio': float(max(r['green_imaginary_ratio'] for r in rows)),
        'max_solver_residual_mV': float(max(r['max_solver_residual_mV'] for r in rows)),
    }

    if aggregate['median_human_to_frozen_high_dose_ratio'] <= 1.05:
        classification = 'NMDA_FEEDBACK_SMALL_IN_REDUCED_ASSAY'
    elif (
        aggregate['human_minus_frozen_effective_rank'] >= 0.5
        and aggregate['human_enhancement_r2_from_Rmed_Rmax_length'] < 0.8
    ):
        classification = 'NMDA_ADDS_BRANCH_DIVERSITY'
    else:
        classification = 'NMDA_ADDS_NONLINEAR_GAIN_NOT_NEW_DIMENSION'

    summary = {
        'gate': 14,
        'object': 'quasi-static NMDA feedback on real reconstructed branch compartments',
        'source_parameters': {
            'human_AMPA_nS': 0.88, 'human_NMDA_nS': 1.31, 'human_gamma_per_mV': 0.078,
            'hybridB_gamma_per_mV': 0.062, 'Mg_mM': 1.0, 'n_per_mM': 1.0/3.57,
        },
        'assay_limits': {
            'quasi_static_peak_conductance': True,
            'does_not_reproduce_synaptic_rise_decay': True,
            'soma_is_voltage_clamped': True,
            'does_not_reproduce_FCI_or_spiking': True,
        },
        'aggregate': aggregate,
        'classification': classification,
        'branches': rows,
        'stopping_line': (
            'This gate asks only whether the released NMDA voltage gate creates '
            'cooperative branch-local responses on the exact passive scaffold. A positive '
            'result must survive the frozen-block and scalar-input-resistance attackers '
            'before a temporal NEURON-level assay is justified.'
        ),
    }
    (out/'gate14.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')

    if not args.no_plot and most_effect is not None:
        _, bi, curves = most_effect
        plot_curves(
            out/'most_nmda_sensitive_branch.png', doses, curves,
            f'Cell 1125 branch {bi} — reduced NMDA feedback',
        )

    print('Operaattori Gate 14 — NMDA closes the branch-local loop')
    print()
    print(f"branches:                              {aggregate['branches']}")
    print(f"median human/frozen at 48 synapses:    {aggregate['median_human_to_frozen_high_dose_ratio']:.4f}")
    print(f"max human/frozen at 48 synapses:       {aggregate['max_human_to_frozen_high_dose_ratio']:.4f}")
    print(f"median human/hybrid-B at 48:            {aggregate['median_human_to_hybridB_high_dose_ratio']:.4f}")
    print(f"branches >10% human/frozen:            {aggregate['fraction_branches_human_over_frozen_10pct']:.3f}")
    print(f"median human cooperativity:            {aggregate['median_human_high_dose_cooperativity']:.4f}")
    print(f"median frozen cooperativity:           {aggregate['median_frozen_high_dose_cooperativity']:.4f}")
    print(f"linear superposition max error:        {aggregate['max_abs_linear_cooperativity_minus_one']:.3e}")
    print(f"dose-curve rank HUMAN:                 {ranks['human']:.3f}")
    print(f"dose-curve rank FROZEN:                {ranks['human_frozen_block']:.3f}")
    print(f"human rank minus frozen:               {aggregate['human_minus_frozen_effective_rank']:.3f}")
    print(f"enhancement R2 from Rmed/Rmax/length:  {scalar_r2:.3f}")
    print(f"max Green symmetry error:              {aggregate['max_green_symmetry_error_Mohm']:.3e} MOhm")
    print(f"max solver residual:                   {aggregate['max_solver_residual_mV']:.3e} mV")
    print()
    print(f'classification: {classification}')

    assert len(rows) >= 6
    assert aggregate['max_green_symmetry_error_Mohm'] < 1e-6
    assert aggregate['max_green_imaginary_ratio'] < 1e-10
    assert aggregate['max_abs_linear_cooperativity_minus_one'] < 1e-9
    assert aggregate['max_solver_residual_mV'] < 1e-5
    assert classification in {
        'NMDA_FEEDBACK_SMALL_IN_REDUCED_ASSAY',
        'NMDA_ADDS_BRANCH_DIVERSITY',
        'NMDA_ADDS_NONLINEAR_GAIN_NOT_NEW_DIMENSION',
    }


if __name__ == '__main__':
    main()
