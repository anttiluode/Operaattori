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

from operaattori.cable_path import (
    PassiveCableParams,
    area_preserving_radius,
    cable_abcd,
    commutator_action_score,
    group_delay_ms,
    path_abcd,
    phase_difference,
    relative_complex_difference,
    sealed_distal_voltage_gain,
    sealed_input_impedance,
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
    req = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Operaattori-Gate11/1.0"})
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
        # The synthetic soma-centroid -> first-neurite edge is geometric
        # scaffold bookkeeping, not a cylindrical dendritic cable segment.
        if parent == 0:
            continue
        length = float(np.linalg.norm(tree.positions[child] - tree.positions[parent]))
        if length <= 1e-8:
            continue
        radius = 0.5 * float(tree.radii[child] + tree.radii[parent])
        radius = max(radius, 0.15)  # paper clamps diameters to >= 0.3 um
        lengths.append(length)
        radii.append(radius)
        child_nodes.append(child)
    return np.asarray(lengths), np.asarray(radii), child_nodes


def plot_sensitive_path(path: Path, freq: np.ndarray, real_h, reverse_h, uniform_h, title: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    ax.semilogx(freq, 20.0 * np.log10(np.maximum(np.abs(real_h), 1e-30)), label="real order")
    ax.semilogx(freq, 20.0 * np.log10(np.maximum(np.abs(reverse_h), 1e-30)), label="reversed same segments")
    ax.semilogx(freq, 20.0 * np.log10(np.maximum(np.abs(uniform_h), 1e-30)), label="area-matched uniform radius")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("sealed distal voltage gain (dB)")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--morphology", type=Path)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results" / "gate11")
    ap.add_argument("--paths", type=int, default=64)
    ap.add_argument("--shuffles", type=int, default=8)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(out / "_source" / SOURCE_NAME)
    tree = load_morphio_tree(source)

    cc = child_counts(tree.parents)
    # MorphIO SectionType.axon == 2; keep dendritic/custom endings and exclude
    # the known incomplete axon from this first cable-path assay.
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
        raise RuntimeError("no usable dendritic tip paths")

    # Dense enough for phase/group-delay measurement, still cheap.
    frequencies = np.geomspace(1.0, 300.0, 36)
    params = PassiveCableParams()
    rng = np.random.default_rng(11025)

    rows = []
    most_sensitive = None
    most_score = -1.0

    for total_len, leaf, lengths, radii, nodes in selected:
        n = len(lengths)
        reverse_order = np.arange(n - 1, -1, -1)
        r_uniform = area_preserving_radius(lengths, radii)
        uniform_radii = np.full_like(radii, r_uniform)

        real_gain = []
        reverse_gain = []
        uniform_gain = []
        real_z = []
        reverse_z = []
        uniform_reverse_z = []
        shuffle_z_diff = []
        adjacent_comm = []

        # A deterministic shuffle family per path, preserving the exact
        # (length,radius) segment multiset.
        perms = [rng.permutation(n) for _ in range(args.shuffles)]

        for f in frequencies:
            Mr = path_abcd(lengths, radii, f, params)
            Mv = path_abcd(lengths, radii, f, params, order=reverse_order)
            Mu = path_abcd(lengths, uniform_radii, f, params)
            Muv = path_abcd(lengths, uniform_radii, f, params, order=reverse_order)

            zr = sealed_input_impedance(Mr)
            zv = sealed_input_impedance(Mv)
            zuv = sealed_input_impedance(Muv)
            real_z.append(zr)
            reverse_z.append(zv)
            uniform_reverse_z.append(relative_complex_difference(sealed_input_impedance(Mu), zuv))

            gr = sealed_distal_voltage_gain(Mr)
            gv = sealed_distal_voltage_gain(Mv)
            gu = sealed_distal_voltage_gain(Mu)
            real_gain.append(gr)
            reverse_gain.append(gv)
            uniform_gain.append(gu)

            zdiffs = []
            for perm in perms:
                Ms = path_abcd(lengths, radii, f, params, order=perm)
                zdiffs.append(
                    relative_complex_difference(zr, sealed_input_impedance(Ms))
                )
            shuffle_z_diff.append(float(np.mean(zdiffs)))

            # Sample adjacent local noncommutativity without making this
            # diagnostic dominate runtime.
            step = max(1, n // 24)
            scores = []
            for j in range(0, n - 1, step):
                A = cable_abcd(lengths[j], radii[j], f, params)
                B = cable_abcd(lengths[j + 1], radii[j + 1], f, params)
                scores.append(commutator_action_score(A, B))
            adjacent_comm.append(float(np.mean(scores)) if scores else 0.0)

        real_gain = np.asarray(real_gain)
        reverse_gain = np.asarray(reverse_gain)
        uniform_gain = np.asarray(uniform_gain)
        real_z = np.asarray(real_z)
        reverse_z = np.asarray(reverse_z)

        z_reverse_diff = np.asarray([
            relative_complex_difference(a, b) for a, b in zip(real_z, reverse_z)
        ])
        gain_reverse_diff = np.asarray([
            relative_complex_difference(a, b) for a, b in zip(real_gain, reverse_gain)
        ])
        phase_reverse_diff = np.asarray([
            phase_difference(a, b) for a, b in zip(real_gain, reverse_gain)
        ])
        gd_real = group_delay_ms(frequencies, real_gain)
        gd_reverse = group_delay_ms(frequencies, reverse_gain)
        gd_diff = np.abs(gd_real - gd_reverse)

        row = {
            "leaf": int(leaf),
            "segments": int(n),
            "path_length_um": float(total_len),
            "radius_min_um": float(np.min(radii)),
            "radius_max_um": float(np.max(radii)),
            "area_preserving_uniform_radius_um": float(r_uniform),
            "median_adjacent_commutator_action": float(np.median(adjacent_comm)),
            "median_reverse_input_impedance_relative_difference": float(np.median(z_reverse_diff)),
            "max_reverse_input_impedance_relative_difference": float(np.max(z_reverse_diff)),
            "median_shuffle_input_impedance_relative_difference": float(np.median(shuffle_z_diff)),
            "median_reverse_gain_relative_difference": float(np.median(gain_reverse_diff)),
            "median_reverse_gain_phase_difference_rad": float(np.median(phase_reverse_diff)),
            "median_reverse_group_delay_difference_ms": float(np.median(gd_diff)),
            "max_uniform_radius_reverse_impedance_difference": float(np.max(uniform_reverse_z)),
        }
        rows.append(row)

        score = row["median_reverse_input_impedance_relative_difference"]
        if score > most_score:
            most_score = score
            most_sensitive = {
                "row": row,
                "real_gain": real_gain,
                "reverse_gain": reverse_gain,
                "uniform_gain": uniform_gain,
            }

    def vals(key: str) -> np.ndarray:
        return np.asarray([r[key] for r in rows], dtype=float)

    summary = {
        "gate": 11,
        "object": "passive transport matrices attached to real scaffold paths",
        "source": {
            "paper": "Aizenbud et al. 2026 PNAS",
            "morphology_identifier": "1125",
            "commit": SOURCE_COMMIT,
            "path": SOURCE_REL,
        },
        "passive_parameters": {
            "Rm_ohm_cm2": params.rm_ohm_cm2,
            "Ra_ohm_cm": params.ra_ohm_cm,
            "Cm_uF_cm2": params.cm_uF_cm2,
        },
        "assay": {
            "dendritic_paths": len(rows),
            "frequencies_hz": frequencies.tolist(),
            "shuffle_permutations_per_path": int(args.shuffles),
            "same_segment_multiset_in_reverse_and_shuffle": True,
            "uniform_control_preserves_path_length_and_surface_area": True,
        },
        "aggregate": {
            "median_path_length_um": float(np.median(vals("path_length_um"))),
            "median_segments_per_path": float(np.median(vals("segments"))),
            "median_adjacent_commutator_action": float(np.median(vals("median_adjacent_commutator_action"))),
            "median_reverse_impedance_difference": float(np.median(vals("median_reverse_input_impedance_relative_difference"))),
            "median_shuffle_impedance_difference": float(np.median(vals("median_shuffle_input_impedance_relative_difference"))),
            "median_reverse_gain_difference": float(np.median(vals("median_reverse_gain_relative_difference"))),
            "median_reverse_gain_phase_difference_rad": float(np.median(vals("median_reverse_gain_phase_difference_rad"))),
            "median_reverse_group_delay_difference_ms": float(np.median(vals("median_reverse_group_delay_difference_ms"))),
            "max_commuting_control_reverse_impedance_difference": float(np.max(vals("max_uniform_radius_reverse_impedance_difference"))),
        },
        "most_order_sensitive_path": most_sensitive["row"] if most_sensitive else None,
        "paths": rows,
        "stopping_line": (
            "Real dendritic geometry supplies a spatially ordered set of passive cable "
            "transport matrices. Radius variation makes those local operators "
            "noncommuting, so reversing or shuffling the exact same segment multiset "
            "changes impedance, attenuation and phase. An area-matched uniform-radius "
            "control commutes and loses this order sensitivity. This earns a physical "
            "path-order degree of freedom, not learning or useful intelligence."
        ),
    }

    (out / "gate11.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    if not args.no_plot and most_sensitive is not None:
        plot_sensitive_path(
            out / "most_order_sensitive_path.png",
            frequencies,
            most_sensitive["real_gain"],
            most_sensitive["reverse_gain"],
            most_sensitive["uniform_gain"],
            (
                f"Cell 1125 tip {most_sensitive['row']['leaf']} — same cable segments, "
                "different spatial order"
            ),
        )

    a = summary["aggregate"]
    print("Operaattori Gate 11 — passive transport phase on the real scaffold")
    print()
    print(f"dendritic paths audited:         {len(rows)}")
    print(f"median path length:             {a['median_path_length_um']:.1f} um")
    print(f"median segments/path:           {a['median_segments_per_path']:.0f}")
    print(f"median adjacent commutator:     {a['median_adjacent_commutator_action']:.3e}")
    print(f"real vs REVERSED impedance:     {a['median_reverse_impedance_difference']:.6g}")
    print(f"real vs SHUFFLED impedance:     {a['median_shuffle_impedance_difference']:.6g}")
    print(f"real vs REVERSED gain:          {a['median_reverse_gain_difference']:.6g}")
    print(f"median phase difference:        {a['median_reverse_gain_phase_difference_rad']:.6g} rad")
    print(f"median group-delay difference:  {a['median_reverse_group_delay_difference_ms']:.6g} ms")
    print(f"COMMUTING control max diff:     {a['max_commuting_control_reverse_impedance_difference']:.3e}")
    print()
    print(summary["stopping_line"])

    # The gate is comparative rather than absolute: the real heterogeneous
    # cable must show measurable order sensitivity, while the matched uniform
    # cable must collapse it to roundoff.
    assert len(rows) >= 16
    assert a["median_adjacent_commutator_action"] > 1e-10
    assert a["median_reverse_impedance_difference"] > 1e-7
    assert a["median_shuffle_impedance_difference"] > 1e-7
    assert a["max_commuting_control_reverse_impedance_difference"] < 1e-9
    assert (
        a["median_reverse_impedance_difference"]
        > 1e3 * a["max_commuting_control_reverse_impedance_difference"]
    )


if __name__ == "__main__":
    main()
