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
from operaattori.order_nulls import (
    OrderAuditConfig,
    compose_orders,
    empirical_unusualness,
    endpoint_preserving_permutations,
    full_permutations,
    linear_taper_r2,
    monotone_radius_orders,
    precompute_segment_matrices,
    standardized_distance,
    transfer_features,
    within_window_permutations,
)
from operaattori.real_scaffold import child_counts, load_morphio_tree


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
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Operaattori-Gate12/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r, dest.open("wb") as f:
        f.write(r.read())
    return dest


def path_nodes_to_root(parents: np.ndarray, leaf: int) -> list[int]:
    out = [int(leaf)]
    i = int(leaf)
    while i > 0:
        i = int(parents[i])
        out.append(i)
    out.reverse()
    return out


def physical_segments(tree, leaf: int) -> tuple[np.ndarray, np.ndarray, list[int]]:
    nodes = path_nodes_to_root(tree.parents, leaf)
    lengths = []
    radii = []
    child_nodes = []
    for child in nodes[1:]:
        parent = int(tree.parents[child])
        if parent == 0:
            continue
        length = float(np.linalg.norm(tree.positions[child] - tree.positions[parent]))
        if length <= 1e-8:
            continue
        radius = 0.5 * float(tree.radii[child] + tree.radii[parent])
        radius = max(radius, 0.15)
        lengths.append(length)
        radii.append(radius)
        child_nodes.append(child)
    return np.asarray(lengths), np.asarray(radii), child_nodes


def scalar_percentile(value: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    return float((1 + np.sum(null <= value)) / (len(null) + 1))


def describe_signal(median_p: float, fraction_p05: float) -> str:
    if median_p <= 0.05 and fraction_p05 >= 0.50:
        return "strong"
    if median_p <= 0.10 or fraction_p05 >= 0.25:
        return "mixed"
    return "weak"


def plot_summary(out: Path, rows: list[dict], windows: tuple[float, ...]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ["full", "endpoint"] + [f"{w:g}um" for w in windows]
    data = [
        [r["full_null"]["empirical_upper_tail_p"] for r in rows],
        [r["endpoint_null"]["empirical_upper_tail_p"] for r in rows],
    ] + [
        [r["coarse_window_nulls"][str(w)]["empirical_upper_tail_p"] for r in rows]
        for w in windows
    ]

    fig = plt.figure(figsize=(9, 5))
    ax = fig.add_subplot(111)
    ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.axhline(0.05, linestyle="--", linewidth=1)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 1.0)
    ax.set_ylabel("empirical upper-tail p (smaller = more unusual)")
    ax.set_xlabel("null model")
    ax.set_title("Cell 1125 — is the biological cable ordering unusual?")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--morphology", type=Path)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results" / "gate12")
    ap.add_argument("--paths", type=int, default=48)
    ap.add_argument("--permutations", type=int, default=256)
    ap.add_argument("--constrained-permutations", type=int, default=96)
    ap.add_argument("--seed", type=int, default=12025)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if args.permutations < 32 or args.constrained_permutations < 24:
        raise ValueError("permutation counts are too small for this audit")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(out / "_source" / SOURCE_NAME)
    tree = load_morphio_tree(source)
    cfg = OrderAuditConfig()
    freqs = np.asarray(cfg.frequencies_hz, dtype=float)
    params = PassiveCableParams()

    cc = child_counts(tree.parents)
    leaves = [
        i for i in range(1, len(tree.parents))
        if cc[i] == 0 and int(tree.section_types[i]) != 2
    ]

    candidates = []
    for leaf in leaves:
        lengths, radii, nodes = physical_segments(tree, leaf)
        if len(lengths) < 20:
            continue
        candidates.append((float(np.sum(lengths)), leaf, lengths, radii, nodes))
    candidates.sort(reverse=True, key=lambda x: x[0])
    selected = candidates[: min(args.paths, len(candidates))]
    if not selected:
        raise RuntimeError("no usable dendritic paths")

    rows = []
    for path_index, (total_len, leaf, lengths, radii, nodes) in enumerate(selected):
        rng = np.random.default_rng(args.seed + 1009 * path_index + int(leaf))
        n = len(lengths)
        mats = precompute_segment_matrices(lengths, radii, freqs, params)

        real_order = np.arange(n, dtype=int)
        real_feat = transfer_features(compose_orders(mats, real_order))[0]

        p_full = full_permutations(n, args.permutations, rng)
        feat_full = transfer_features(compose_orders(mats, p_full))
        full_stats = empirical_unusualness(real_feat, feat_full)

        thick_to_thin, thin_to_thick = monotone_radius_orders(radii)
        feat_mono = transfer_features(
            compose_orders(mats, np.stack([thick_to_thin, thin_to_thick]))
        )
        thick_stats = empirical_unusualness(feat_mono[0], feat_full)
        thin_stats = empirical_unusualness(feat_mono[1], feat_full)

        p_endpoint = endpoint_preserving_permutations(
            lengths,
            args.constrained_permutations,
            rng,
            fraction=cfg.endpoint_fraction,
        )
        feat_endpoint = transfer_features(compose_orders(mats, p_endpoint))
        endpoint_stats = empirical_unusualness(real_feat, feat_endpoint)

        coarse = {}
        for window in cfg.coarse_windows_um:
            p_local = within_window_permutations(
                lengths,
                window,
                args.constrained_permutations,
                rng,
            )
            feat_local = transfer_features(compose_orders(mats, p_local))
            coarse[str(window)] = empirical_unusualness(real_feat, feat_local)

        slope, taper_r2 = linear_taper_r2(lengths, radii)

        F = len(freqs)
        real_gain_db = real_feat[:F]
        null_gain_db = feat_full[:, :F]
        gain_percentiles = {
            str(float(f)): scalar_percentile(float(real_gain_db[j]), null_gain_db[:, j])
            for j, f in enumerate(freqs)
        }

        radius_increases = np.diff(radii) > 1e-12
        thick_is_real = bool(np.array_equal(thick_to_thin, real_order))
        real_thick_max_abs_feature_error = float(np.max(np.abs(real_feat - feat_mono[0])))

        row = {
            "leaf": int(leaf),
            "segments": int(n),
            "path_length_um": float(total_len),
            "radius_min_um": float(np.min(radii)),
            "radius_max_um": float(np.max(radii)),
            "linear_taper_slope_um_per_normalized_path": float(slope),
            "linear_taper_r2": float(taper_r2),
            "radius_increase_steps": int(np.sum(radius_increases)),
            "fraction_radius_steps_nonincreasing": float(np.mean(~radius_increases)) if len(radius_increases) else 1.0,
            "stable_thick_to_thin_order_is_exact_real_order": thick_is_real,
            "real_vs_thick_to_thin_max_abs_transfer_feature_error": real_thick_max_abs_feature_error,
            "full_null": full_stats,
            "endpoint_null": endpoint_stats,
            "coarse_window_nulls": coarse,
            "thick_to_thin_null_position": thick_stats,
            "thin_to_thick_null_position": thin_stats,
            "real_to_thick_to_thin_standardized_distance": standardized_distance(
                real_feat, feat_mono[0], feat_full
            ),
            "real_to_thin_to_thick_standardized_distance": standardized_distance(
                real_feat, feat_mono[1], feat_full
            ),
            "real_gain_db_percentiles_under_full_null": gain_percentiles,
        }
        rows.append(row)

        print(
            f"[{path_index+1:02d}/{len(selected):02d}] leaf {leaf:5d} "
            f"L={total_len:7.1f}um n={n:4d} "
            f"full p={full_stats['empirical_upper_tail_p']:.4f} "
            f"endpoint p={endpoint_stats['empirical_upper_tail_p']:.4f}"
        )

    def arr(fn):
        return np.asarray([fn(r) for r in rows], dtype=float)

    full_p = arr(lambda r: r["full_null"]["empirical_upper_tail_p"])
    endpoint_p = arr(lambda r: r["endpoint_null"]["empirical_upper_tail_p"])

    coarse_agg = {}
    for w in cfg.coarse_windows_um:
        p = arr(lambda r, w=w: r["coarse_window_nulls"][str(w)]["empirical_upper_tail_p"])
        coarse_agg[str(w)] = {
            "median_tail_p": float(np.median(p)),
            "fraction_tail_p_le_0.05": float(np.mean(p <= 0.05)),
            "signal": describe_signal(float(np.median(p)), float(np.mean(p <= 0.05))),
        }

    aggregate = {
        "median_path_length_um": float(np.median(arr(lambda r: r["path_length_um"]))),
        "median_segments_per_path": float(np.median(arr(lambda r: r["segments"]))),
        "median_linear_taper_r2": float(np.median(arr(lambda r: r["linear_taper_r2"]))),
        "median_linear_taper_slope": float(np.median(arr(lambda r: r["linear_taper_slope_um_per_normalized_path"]))),
        "median_radius_increase_steps": float(np.median(arr(lambda r: r["radius_increase_steps"]))),
        "median_fraction_radius_steps_nonincreasing": float(np.median(arr(lambda r: r["fraction_radius_steps_nonincreasing"]))),
        "fraction_paths_exactly_equal_stable_thick_to_thin": float(np.mean(arr(lambda r: float(r["stable_thick_to_thin_order_is_exact_real_order"])))),
        "max_real_vs_thick_to_thin_transfer_feature_error": float(np.max(arr(lambda r: r["real_vs_thick_to_thin_max_abs_transfer_feature_error"]))),
        "full_permutation": {
            "median_tail_p": float(np.median(full_p)),
            "fraction_tail_p_le_0.05": float(np.mean(full_p <= 0.05)),
            "signal": describe_signal(float(np.median(full_p)), float(np.mean(full_p <= 0.05))),
        },
        "endpoint_preserving": {
            "median_tail_p": float(np.median(endpoint_p)),
            "fraction_tail_p_le_0.05": float(np.mean(endpoint_p <= 0.05)),
            "signal": describe_signal(float(np.median(endpoint_p)), float(np.mean(endpoint_p <= 0.05))),
        },
        "coarse_window_preserving": coarse_agg,
        "median_real_to_thick_to_thin_distance_z": float(np.median(arr(lambda r: r["real_to_thick_to_thin_standardized_distance"]))),
        "median_real_to_thin_to_thick_distance_z": float(np.median(arr(lambda r: r["real_to_thin_to_thick_standardized_distance"]))),
        "median_thick_to_thin_tail_p": float(np.median(arr(lambda r: r["thick_to_thin_null_position"]["empirical_upper_tail_p"]))),
        "median_thin_to_thick_tail_p": float(np.median(arr(lambda r: r["thin_to_thick_null_position"]["empirical_upper_tail_p"]))),
    }

    full_signal = aggregate["full_permutation"]["signal"]
    endpoint_signal = aggregate["endpoint_preserving"]["signal"]
    p50_signal = aggregate["coarse_window_preserving"]["50.0"]["signal"]

    if (
        aggregate["fraction_paths_exactly_equal_stable_thick_to_thin"] >= 0.90
        or aggregate["median_real_to_thick_to_thin_distance_z"] <= 1e-6
    ):
        classification = "MONOTONIC_TAPER_EXPLAINS_REAL_ORDER"
        interpretation = (
            "The biological paths are already the same ordering produced by the simple "
            "stable thick-to-thin taper ruler (or numerically indistinguishable from it). "
            "The permutation tail result therefore does not establish additional fine "
            "operator ordering beyond ordinary taper."
        )
    elif full_signal == "weak":
        classification = "REAL_ORDER_NOT_GLOBALLY_UNUSUAL"
        interpretation = (
            "Gate 11's order sensitivity is real, but the biological sequence is not "
            "consistently in the tail of the unconstrained permutation null."
        )
    elif endpoint_signal == "weak":
        classification = "ENDPOINT_POSITION_EXPLAINS_MUCH"
        interpretation = (
            "The biological transfer is unusual under full permutations but becomes "
            "ordinary when proximal/distal zones are preserved."
        )
    elif p50_signal == "weak":
        classification = "COARSE_TAPER_EXPLAINS_MUCH"
        interpretation = (
            "The biological transfer is unusual under broad reorderings but becomes "
            "ordinary when the coarse 50-um proximal-distal profile is preserved."
        )
    else:
        classification = "FINE_ORDER_REMAINS_UNUSUAL"
        interpretation = (
            "The biological transfer remains unusually positioned even against "
            "endpoint- and coarse-profile-preserving nulls. This is an order-tuning "
            "signal, not evidence of optimality."
        )

    summary = {
        "gate": 12,
        "object": "biological cable-order audit on real human L2/3 cell 1125",
        "source": {
            "paper": "Aizenbud et al. 2026 PNAS",
            "morphology_identifier": "1125",
            "commit": SOURCE_COMMIT,
            "path": SOURCE_REL,
        },
        "assay": {
            "paths": len(rows),
            "full_permutations_per_path": int(args.permutations),
            "constrained_permutations_per_path": int(args.constrained_permutations),
            "frequencies_hz": freqs.tolist(),
            "endpoint_fraction_each_side": cfg.endpoint_fraction,
            "coarse_windows_um": list(cfg.coarse_windows_um),
            "path_level_results_are_not_independent_cell_replicates": True,
            "no_task_utility_or_optimality_claim": True,
        },
        "aggregate": aggregate,
        "classification": classification,
        "interpretation": interpretation,
        "paths": rows,
        "stopping_line": (
            "Gate 12 distinguishes mere order sensitivity from biological-order "
            "unusualness. A positive tail result still does not mean the ordering is "
            "better or evolved for the measured passive transfer."
        ),
    }

    (out / "gate12.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if not args.no_plot:
        plot_summary(out / "biological_order_nulls.png", rows, cfg.coarse_windows_um)

    print()
    print("Operaattori Gate 12 — is the biological cable order special?")
    print()
    print(f"paths:                           {len(rows)}")
    print(f"full permutations/path:          {args.permutations}")
    print(f"constrained permutations/path:   {args.constrained_permutations}")
    print(f"median gross taper R^2:          {aggregate['median_linear_taper_r2']:.4f}")
    print(f"median radius-increase steps:    {aggregate['median_radius_increase_steps']:.1f}")
    print(f"median nonincreasing fraction:   {aggregate['median_fraction_radius_steps_nonincreasing']:.6f}")
    print(f"paths exactly thick->thin:       {aggregate['fraction_paths_exactly_equal_stable_thick_to_thin']:.3f}")
    print(f"max real/thick transfer error:   {aggregate['max_real_vs_thick_to_thin_transfer_feature_error']:.3e}")
    print(
        "FULL null median p / p<=.05:    "
        f"{aggregate['full_permutation']['median_tail_p']:.4f} / "
        f"{aggregate['full_permutation']['fraction_tail_p_le_0.05']:.3f}"
    )
    print(
        "ENDPOINT null median p / p<=.05:"
        f" {aggregate['endpoint_preserving']['median_tail_p']:.4f} / "
        f"{aggregate['endpoint_preserving']['fraction_tail_p_le_0.05']:.3f}"
    )
    for w in cfg.coarse_windows_um:
        a = aggregate["coarse_window_preserving"][str(w)]
        print(
            f"{w:>5g}um null median p / p<=.05: "
            f"{a['median_tail_p']:.4f} / {a['fraction_tail_p_le_0.05']:.3f}"
        )
    print(
        "real -> thick-to-thin distance: "
        f"{aggregate['median_real_to_thick_to_thin_distance_z']:.3f} null-SD RMS"
    )
    print(
        "real -> thin-to-thick distance: "
        f"{aggregate['median_real_to_thin_to_thick_distance_z']:.3f} null-SD RMS"
    )
    print()
    print(f"classification: {classification}")
    print(interpretation)

    assert len(rows) >= 16
    assert np.all(np.isfinite(full_p))
    assert np.all((full_p > 0.0) & (full_p <= 1.0))
    assert aggregate["median_segments_per_path"] >= 20
    assert classification in {
        "MONOTONIC_TAPER_EXPLAINS_REAL_ORDER",
        "REAL_ORDER_NOT_GLOBALLY_UNUSUAL",
        "ENDPOINT_POSITION_EXPLAINS_MUCH",
        "COARSE_TAPER_EXPLAINS_MUCH",
        "FINE_ORDER_REMAINS_UNUSUAL",
    }


if __name__ == "__main__":
    main()
