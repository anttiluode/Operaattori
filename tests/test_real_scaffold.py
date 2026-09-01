import unittest

import numpy as np

from operaattori.real_scaffold import (
    PointTree,
    build_matrix_scaffold,
    descendant_mask,
    edge_lengths,
    reconstruct,
    rotation_quality,
    twist_scaffold,
)


def synthetic_tree() -> PointTree:
    # root -> trunk -> bifurcation -> two bent daughters
    p = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 2.0],
            [1.0, 0.0, 3.0],
            [2.0, 0.5, 4.0],
            [-1.0, 0.0, 3.0],
            [-2.0, -0.5, 4.0],
        ]
    )
    parents = np.asarray([-1, 0, 1, 2, 3, 2, 5], dtype=np.int64)
    n = len(p)
    return PointTree(
        positions=p,
        parents=parents,
        radii=np.ones(n),
        section_ids=np.arange(n, dtype=np.int64),
        section_types=np.full(n, 3, dtype=np.int64),
        soma_points=np.asarray([[0.0, 0.0, 0.0]]),
        soma_radii=np.asarray([1.0]),
    )


class RealScaffoldTests(unittest.TestCase):
    def test_local_matrices_reconstruct_absolute_geometry(self):
        tree = synthetic_tree()
        scaffold = build_matrix_scaffold(tree)
        recon, _ = reconstruct(scaffold)
        self.assertLess(float(np.max(np.abs(recon - tree.positions))), 1e-12)

    def test_local_rotations_are_rigid_frames(self):
        scaffold = build_matrix_scaffold(synthetic_tree())
        q = rotation_quality(scaffold.local_transforms)
        self.assertLess(q["max_orthogonality_error"], 1e-12)
        self.assertLess(q["max_abs_det_minus_one"], 1e-12)

    def test_one_local_bend_moves_only_distal_subtree(self):
        tree = synthetic_tree()
        scaffold = build_matrix_scaffold(tree)
        p0, _ = reconstruct(scaffold)
        bent = twist_scaffold(scaffold, pivot=3, angle_degrees=25.0, axis="y")
        p1, _ = reconstruct(bent)

        mask = descendant_mask(tree.parents, 3)
        disp = np.linalg.norm(p1 - p0, axis=1)

        self.assertLess(disp[3], 1e-12)  # attachment point itself stays fixed
        self.assertGreater(float(np.max(disp[mask])), 0.1)
        self.assertLess(float(np.max(disp[~mask])), 1e-12)

        l0 = edge_lengths(p0, tree.parents)
        l1 = edge_lengths(p1, tree.parents)
        self.assertLess(float(np.max(np.abs(l0 - l1))), 1e-12)


if __name__ == "__main__":
    unittest.main()
