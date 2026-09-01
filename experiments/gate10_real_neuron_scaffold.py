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

from operaattori.real_scaffold import (
    build_matrix_scaffold,
    choose_twist_pivot,
    descendant_counts,
    descendant_mask,
    edge_lengths,
    export_npz,
    load_morphio_tree,
    plot_tree,
    reconstruct,
    rotation_quality,
    twist_scaffold,
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
    print(f"Downloading pinned Aizenbud morphology:\n  {SOURCE_URL}")
    req = urllib.request.Request(
        SOURCE_URL,
        headers={"User-Agent": "Operaattori-Gate10/1.0"},
    )
    with urllib.request.urlopen(req, timeout=60) as r, dest.open("wb") as f:
        f.write(r.read())
    return dest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--morphology",
        type=Path,
        help="optional local ASC/SWC/H5; otherwise download the pinned Aizenbud cell",
    )
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results" / "gate10")
    ap.add_argument("--twist-deg", type=float, default=20.0)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    source = args.morphology or download_source(out / "_source" / SOURCE_NAME)

    tree = load_morphio_tree(source)
    scaffold = build_matrix_scaffold(tree)
    recon, _ = reconstruct(scaffold)

    err = np.linalg.norm(recon - tree.positions, axis=1)
    original_lengths = edge_lengths(tree.positions, tree.parents)
    recon_lengths = edge_lengths(recon, tree.parents)

    pivot = choose_twist_pivot(tree.parents)
    twisted = twist_scaffold(scaffold, pivot, angle_degrees=args.twist_deg, axis="y")
    twisted_pos, _ = reconstruct(twisted)
    mask = descendant_mask(tree.parents, pivot)

    disp = np.linalg.norm(twisted_pos - recon, axis=1)
    outside = ~mask
    twisted_lengths = edge_lengths(twisted_pos, tree.parents)
    dcounts = descendant_counts(tree.parents)

    child_count = np.zeros(len(tree.parents), dtype=np.int64)
    for i in range(1, len(tree.parents)):
        child_count[int(tree.parents[i])] += 1

    rq = rotation_quality(scaffold.local_transforms)

    summary = {
        "gate": 10,
        "object": "branching local SE(3) scaffold on a real human cortical neuron",
        "source": {
            "paper": "Aizenbud et al. 2026 PNAS",
            "repository": "ido4848/FCI",
            "commit": SOURCE_COMMIT,
            "path": SOURCE_REL,
            "morphology_identifier": "1125",
            "cell": "human L2/3 pyramidal neuron",
            "reported_fci": 0.4294,
        },
        "morphology": {
            "nodes_including_synthetic_soma_root": int(len(tree.positions)),
            "edges": int(len(tree.positions) - 1),
            "soma_contour_points": int(len(tree.soma_points)),
            "branch_nodes": int(np.sum(child_count >= 2)),
            "tips": int(np.sum(child_count == 0)),
            "total_cable_length_um": float(np.sum(original_lengths)),
            "bbox_um": (np.max(tree.positions, axis=0) - np.min(tree.positions, axis=0)).tolist(),
            "section_type_counts": {
                str(int(k)): int(v)
                for k, v in zip(*np.unique(tree.section_types, return_counts=True))
            },
        },
        "matrix_scaffold": {
            "stored_absolute_neurite_coordinates": False,
            "root_absolute_coordinates": 3,
            "local_matrix_shape": [4, 4],
            "local_matrices": int(len(scaffold.local_transforms) - 1),
            "max_reconstruction_error_um": float(np.max(err)),
            "rms_reconstruction_error_um": float(np.sqrt(np.mean(err * err))),
            "max_edge_length_error_um": float(np.max(np.abs(recon_lengths - original_lengths))),
            **rq,
        },
        "single_local_twist": {
            "pivot_node": int(pivot),
            "pivot_descendants": int(dcounts[pivot]),
            "angle_degrees": float(args.twist_deg),
            "pivot_displacement_um": float(disp[pivot]),
            "max_distal_displacement_um": float(np.max(disp[mask])),
            "mean_distal_displacement_um": float(np.mean(disp[mask])),
            "max_outside_subtree_displacement_um": float(np.max(disp[outside])) if np.any(outside) else 0.0,
            "max_edge_length_change_um": float(np.max(np.abs(twisted_lengths - recon_lengths))),
        },
        "stopping_line": (
            "The real arbor can be represented exactly as a branching field of parent-local "
            "SE(3) matrices. A change to one local frame coherently moves only its distal "
            "subtree while preserving cable lengths. This earns geometry-as-scaffold, not "
            "learning or neuronal superiority. Gate 11 may attach dynamical operators to "
            "these local frames and test whether spatial operator order matters."
        ),
    }

    (out / "gate10.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    export_npz(out / "cell1125_matrix_scaffold.npz", scaffold)

    if not args.no_plot:
        plot_tree(
            out / "cell1125_original_scaffold.png",
            recon,
            tree.parents,
            title="Aizenbud human L2/3 cell 1125 — reconstructed from local matrices",
        )
        plot_tree(
            out / "cell1125_one_local_bend.png",
            twisted_pos,
            tree.parents,
            title=f"Same scaffold after one local {args.twist_deg:g} degree bend",
            pivot=pivot,
        )

    print("Operaattori Gate 10 — real neuron matrix scaffold")
    print()
    print(f"source:                     {source}")
    print(f"nodes / edges:              {len(tree.positions)} / {len(tree.positions)-1}")
    print(f"branch nodes / tips:        {summary['morphology']['branch_nodes']} / {summary['morphology']['tips']}")
    print(f"total cable length:         {summary['morphology']['total_cable_length_um']:.2f} um")
    print(f"max reconstruction error:  {summary['matrix_scaffold']['max_reconstruction_error_um']:.3e} um")
    print(f"max rotation orth error:    {rq['max_orthogonality_error']:.3e}")
    print()
    print(f"twist pivot:                {pivot}")
    print(f"distal descendants:         {dcounts[pivot]}")
    print(f"pivot displacement:         {disp[pivot]:.3e} um")
    print(f"max distal displacement:    {np.max(disp[mask]):.3f} um")
    print(f"max outside displacement:   {np.max(disp[outside]) if np.any(outside) else 0.0:.3e} um")
    print(f"max cable-length change:    {np.max(np.abs(twisted_lengths-recon_lengths)):.3e} um")
    print()
    print(summary["stopping_line"])

    # Gate conditions: representation must be nontrivial and exact; a local
    # rotation must bend a sizeable real subtree without teleporting the pivot
    # or changing cable lengths.
    assert len(tree.positions) > 1000
    assert summary["morphology"]["branch_nodes"] > 20
    assert summary["matrix_scaffold"]["max_reconstruction_error_um"] < 1e-7
    assert rq["max_orthogonality_error"] < 1e-10
    assert rq["max_abs_det_minus_one"] < 1e-10
    assert dcounts[pivot] >= 20
    assert disp[pivot] < 1e-8
    assert np.max(disp[mask]) > 1.0
    assert (np.max(disp[outside]) if np.any(outside) else 0.0) < 1e-7
    assert np.max(np.abs(twisted_lengths - recon_lengths)) < 1e-7


if __name__ == "__main__":
    main()
