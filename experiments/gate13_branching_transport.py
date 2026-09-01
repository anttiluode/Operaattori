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

from operaattori.cable_path import PassiveCableParams
from operaattori.real_scaffold import load_morphio_tree
from operaattori.tree_cable import (
    active_child_counts,
    isolated_path_transfer,
    path_to_clamp,
    path_with_side_shunts_transfer,
    side_shunt_admittances,
    solve_tree_frequency,
    standardized_effective_rank,
    transfer_signature_features,
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
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Operaattori-Gate13/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r, dest.open("wb") as f:
        f.write(r.read())
    return dest


def build_electrical_tree(tree):
    """Collapse the synthetic soma-centroid links into voltage-clamped roots."""
    n = len(tree.parents)
    parents = tree.parents.copy()

    # Keep soma bookkeeping node plus basal/apical dendrites. Exclude the known
    # incomplete axon (MorphIO/Neurolucida type 2).
    active = np.zeros(n, dtype=bool)
    active[0] = True
    active[1:] = tree.section_types[1:] != 2

    clamped = np.zeros(n, dtype=bool)
    clamped[0] = True
    for i in range(1, n):
        if active[i] and int(parents[i]) == 0:
            # First point of each dendritic root section is an ideal somatic
            # voltage-clamp boundary. The synthetic soma-centroid edge is not
            # treated as cable.
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


def wrapped_phase_difference(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.abs(np.angle(a * np.conj(b)))


def relative_complex_difference(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.abs(a - b) / (0.5 * (np.abs(a) + np.abs(b)) + 1e-300)


def plot_tip(out: Path, frequencies, full, isolated, tip: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    ax.semilogx(
        frequencies,
        20.0 * np.log10(np.maximum(np.abs(full), 1e-300)),
        label="full branching tree",
    )
    ax.semilogx(
        frequencies,
        20.0 * np.log10(np.maximum(np.abs(isolated), 1e-300)),
        label="same soma-tip path, side branches removed",
    )
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("|I_soma / I_injected| (dB)")
    ax.set_title(f"Cell 1125 tip {tip} — branch loading")
    ax.legend()
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--morphology", type=Path)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results" / "gate13")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(out / "_source" / SOURCE_NAME)
    tree = load_morphio_tree(source)

    parents, lengths, radii, active, clamped = build_electrical_tree(tree)
    child_count = active_child_counts(parents, active)
    tips = np.asarray(
        [
            i for i in range(1, len(parents))
            if active[i] and not clamped[i] and child_count[i] == 0
        ],
        dtype=int,
    )
    if len(tips) < 20:
        raise RuntimeError("too few dendritic tips after electrical-tree filtering")

    frequencies = np.asarray([1.0, 5.0, 15.0, 40.0, 100.0, 300.0])
    params = PassiveCableParams()

    full = np.zeros((len(tips), len(frequencies)), dtype=np.complex128)
    isolated = np.zeros_like(full)
    via_shunts = np.zeros_like(full)
    branch_counts = np.zeros(len(tips), dtype=int)
    path_lengths = np.zeros(len(tips), dtype=float)

    for fi, freq in enumerate(frequencies):
        state = solve_tree_frequency(
            parents, lengths, radii, active, clamped, float(freq), params
        )

        for ti, tip in enumerate(tips):
            path = path_to_clamp(parents, clamped, int(tip))
            if fi == 0:
                path_lengths[ti] = float(np.sum(lengths[path]))
                branch_counts[ti] = int(
                    np.sum([child_count[node] >= 2 for node in path])
                )

            full[ti, fi] = state.transfer_to_clamp[tip]
            isolated[ti, fi] = isolated_path_transfer(
                path, lengths, radii, float(freq), params
            )

            shunts = side_shunt_admittances(
                parents, active, clamped, path, state
            )
            via_shunts[ti, fi] = path_with_side_shunts_transfer(
                path, lengths, radii, shunts, float(freq), params
            )

    representation_rel_error = relative_complex_difference(full, via_shunts)
    branch_effect = relative_complex_difference(full, isolated)
    gain_delta_db = (
        20.0 * np.log10(np.maximum(np.abs(full), 1e-300))
        - 20.0 * np.log10(np.maximum(np.abs(isolated), 1e-300))
    )
    phase_delta = wrapped_phase_difference(full, isolated)

    per_tip_effect = np.median(branch_effect, axis=1)
    if np.std(branch_counts) > 0 and np.std(per_tip_effect) > 0:
        branch_count_corr = float(np.corrcoef(branch_counts, per_tip_effect)[0, 1])
    else:
        branch_count_corr = 0.0

    full_features = transfer_signature_features(full)
    path_features = transfer_signature_features(isolated)
    full_rank = standardized_effective_rank(full_features)
    path_rank = standardized_effective_rank(path_features)

    most_idx = int(np.argmax(per_tip_effect))
    most_tip = int(tips[most_idx])

    aggregate = {
        "dendritic_tips": int(len(tips)),
        "median_path_length_um": float(np.median(path_lengths)),
        "median_branch_junctions_per_tip_path": float(np.median(branch_counts)),
        "max_tree_vs_shunt_product_relative_error": float(np.max(representation_rel_error)),
        "median_full_vs_isolated_relative_difference": float(np.median(branch_effect)),
        "median_abs_gain_change_db": float(np.median(np.abs(gain_delta_db))),
        "median_signed_gain_change_db": float(np.median(gain_delta_db)),
        "median_phase_change_rad": float(np.median(phase_delta)),
        "fraction_tip_frequency_points_over_10pct_change": float(np.mean(branch_effect > 0.10)),
        "fraction_tips_median_effect_over_10pct": float(np.mean(per_tip_effect > 0.10)),
        "branch_count_vs_effect_correlation": branch_count_corr,
        "full_tree_signature_effective_rank": float(full_rank),
        "isolated_path_signature_effective_rank": float(path_rank),
        "effective_rank_difference": float(full_rank - path_rank),
    }

    if aggregate["median_full_vs_isolated_relative_difference"] < 1e-3:
        classification = "BRANCH_LOADING_NEGLIGIBLE"
    elif abs(aggregate["effective_rank_difference"]) >= 0.5:
        classification = "BRANCH_LOADING_CHANGES_TRANSFER_AND_PORTFOLIO"
    else:
        classification = "BRANCH_LOADING_CHANGES_TRANSFER"

    rows = []
    for ti, tip in enumerate(tips):
        rows.append(
            {
                "tip": int(tip),
                "path_length_um": float(path_lengths[ti]),
                "branch_junctions_on_path": int(branch_counts[ti]),
                "median_full_vs_isolated_relative_difference": float(per_tip_effect[ti]),
                "median_abs_gain_change_db": float(np.median(np.abs(gain_delta_db[ti]))),
                "median_phase_change_rad": float(np.median(phase_delta[ti])),
            }
        )

    summary = {
        "gate": 13,
        "object": "full branching passive-cable operator scaffold",
        "source": {
            "paper": "Aizenbud et al. 2026 PNAS",
            "morphology_identifier": "1125",
            "commit": SOURCE_COMMIT,
            "path": SOURCE_REL,
        },
        "assay": {
            "frequencies_hz": frequencies.tolist(),
            "axon_excluded": True,
            "synthetic_soma_centroid_edges_excluded": True,
            "root_dendritic_section_points_voltage_clamped": True,
            "side_branch_ablation_preserves_the_same_soma_tip_path": True,
        },
        "aggregate": aggregate,
        "classification": classification,
        "most_branch_loaded_tip": rows[most_idx],
        "tips": rows,
        "stopping_line": (
            "Gate 13 asks only whether branch junctions add function beyond the same "
            "serial tapering paths. Whole-tree elimination must agree with an explicit "
            "cable-plus-side-shunt matrix product before interpreting the branch effect."
        ),
    }

    (out / "gate13.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if not args.no_plot:
        plot_tip(
            out / "most_branch_loaded_tip.png",
            frequencies,
            full[most_idx],
            isolated[most_idx],
            most_tip,
        )

    print("Operaattori Gate 13 — branching matrices on the real neuron")
    print()
    print(f"dendritic tips:                    {aggregate['dendritic_tips']}")
    print(f"median path length:                {aggregate['median_path_length_um']:.1f} um")
    print(f"median branch junctions/path:      {aggregate['median_branch_junctions_per_tip_path']:.1f}")
    print(f"tree vs shunt-product max error:   {aggregate['max_tree_vs_shunt_product_relative_error']:.3e}")
    print(f"full vs isolated median diff:      {aggregate['median_full_vs_isolated_relative_difference']:.4f}")
    print(f"median |gain change|:              {aggregate['median_abs_gain_change_db']:.3f} dB")
    print(f"median signed gain change:         {aggregate['median_signed_gain_change_db']:.3f} dB")
    print(f"median phase change:               {aggregate['median_phase_change_rad']:.4f} rad")
    print(f"tip-frequency >10% changed:        {aggregate['fraction_tip_frequency_points_over_10pct_change']:.3f}")
    print(f"tips median >10% changed:          {aggregate['fraction_tips_median_effect_over_10pct']:.3f}")
    print(f"branch-count/effect corr:          {aggregate['branch_count_vs_effect_correlation']:.3f}")
    print(f"signature effective rank FULL:     {full_rank:.3f}")
    print(f"signature effective rank PATH:     {path_rank:.3f}")
    print()
    print(f"classification: {classification}")

    # Numerical representation gate, not a biological-significance threshold.
    assert len(tips) >= 20
    assert aggregate["max_tree_vs_shunt_product_relative_error"] < 1e-9
    assert np.all(np.isfinite(branch_effect))
    assert classification in {
        "BRANCH_LOADING_NEGLIGIBLE",
        "BRANCH_LOADING_CHANGES_TRANSFER",
        "BRANCH_LOADING_CHANGES_TRANSFER_AND_PORTFOLIO",
    }


if __name__ == "__main__":
    main()
