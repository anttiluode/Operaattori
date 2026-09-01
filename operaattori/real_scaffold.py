"""Real-neuron matrix scaffold for Operaattori Gate 10.

The scaffold is intentionally geometric before it is adaptive.

A rooted morphology is converted into local SE(3) transforms.  Absolute
coordinates are then discarded; the neuron is reconstructed only by composing
the local transforms from the soma outward.

For node i with parent p:

    W_i = W_p L_i

where W is the world pose and L_i is the parent-local transform.  At a
bifurcation the same parent pose feeds multiple independent child transforms,
so the scaffold literally splits with the arbor.

Changing only the rotation block of one local transform leaves the attachment
point fixed but coherently moves the complete distal subtree.  This is the
minimal "matrix scaffold living on the neuron" object needed before adding
signal dynamics or growth.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class PointTree:
    positions: np.ndarray
    parents: np.ndarray
    radii: np.ndarray
    section_ids: np.ndarray
    section_types: np.ndarray
    soma_points: np.ndarray
    soma_radii: np.ndarray

    def validate(self) -> None:
        n = len(self.positions)
        if self.positions.shape != (n, 3):
            raise ValueError("positions must be N x 3")
        for a in (self.parents, self.radii, self.section_ids, self.section_types):
            if len(a) != n:
                raise ValueError("point-tree arrays must have equal length")
        if n == 0 or int(self.parents[0]) != -1:
            raise ValueError("node 0 must be the synthetic soma root")
        for i in range(1, n):
            p = int(self.parents[i])
            if not 0 <= p < i:
                raise ValueError("parents must precede children")


@dataclass(frozen=True)
class MatrixScaffold:
    root_position: np.ndarray
    parents: np.ndarray
    local_transforms: np.ndarray
    radii: np.ndarray
    section_ids: np.ndarray
    section_types: np.ndarray

    def validate(self) -> None:
        n = len(self.parents)
        if self.local_transforms.shape != (n, 4, 4):
            raise ValueError("local transforms must be N x 4 x 4")
        if self.root_position.shape != (3,):
            raise ValueError("root_position must have shape (3,)")
        if int(self.parents[0]) != -1:
            raise ValueError("root parent must be -1")


def _enum_int(value: object) -> int:
    try:
        return int(value)  # pybind enums normally support this
    except (TypeError, ValueError):
        raw = getattr(value, "value", None)
        if raw is None:
            return -1
        return int(raw)


def load_morphio_tree(path: str | Path, duplicate_tol: float = 1e-5) -> PointTree:
    """Load ASC/SWC/H5 through MorphIO and flatten sections into a point tree.

    MorphIO owns the file-format complexity.  This function only converts its
    section tree into an explicit parent-before-child point tree.

    The soma is represented by one synthetic root at its centroid.  Original
    soma contour points are retained separately for plotting/inspection.
    """
    from morphio import Morphology

    morph = Morphology(str(path))

    soma_points = np.asarray(morph.soma.points, dtype=float)
    soma_diameters = np.asarray(morph.soma.diameters, dtype=float)
    if soma_points.ndim != 2 or (soma_points.size and soma_points.shape[1] != 3):
        raise ValueError("unexpected soma point array")

    if len(soma_points):
        root = np.mean(soma_points, axis=0)
        root_radius = float(np.mean(soma_diameters) / 2.0) if len(soma_diameters) else 1.0
    else:
        first = [np.asarray(s.points[0], dtype=float) for s in morph.root_sections if len(s.points)]
        if not first:
            raise ValueError("morphology has neither soma nor neurite points")
        root = np.mean(np.stack(first), axis=0)
        root_radius = 1.0

    positions: list[np.ndarray] = [root]
    parents: list[int] = [-1]
    radii: list[float] = [root_radius]
    section_ids: list[int] = [-1]
    section_types: list[int] = [1]  # conventional soma code; metadata only

    section_end: dict[int, int] = {}

    # MorphIO's default morphology iterator is depth-first and parent-first.
    for section in morph.iter():
        sid = int(section.id)
        stype = _enum_int(section.type)
        if section.is_root:
            current_parent = 0
        else:
            pid = int(section.parent.id)
            if pid not in section_end:
                raise ValueError(f"parent section {pid} was not visited before child {sid}")
            current_parent = section_end[pid]

        pts = np.asarray(section.points, dtype=float)
        diams = np.asarray(section.diameters, dtype=float)
        if len(pts) != len(diams):
            raise ValueError(f"section {sid}: point/diameter length mismatch")
        if not len(pts):
            section_end[sid] = current_parent
            continue

        start = 0
        if np.linalg.norm(pts[0] - positions[current_parent]) <= duplicate_tol:
            # Morphology formats commonly repeat the parent bifurcation point
            # as the first point of every daughter section.
            start = 1

        for j in range(start, len(pts)):
            p = np.asarray(pts[j], dtype=float)
            if np.linalg.norm(p - positions[current_parent]) <= duplicate_tol:
                continue
            idx = len(positions)
            positions.append(p)
            parents.append(current_parent)
            radii.append(float(diams[j]) / 2.0)
            section_ids.append(sid)
            section_types.append(stype)
            current_parent = idx

        section_end[sid] = current_parent

    tree = PointTree(
        positions=np.asarray(positions, dtype=float),
        parents=np.asarray(parents, dtype=np.int64),
        radii=np.asarray(radii, dtype=float),
        section_ids=np.asarray(section_ids, dtype=np.int64),
        section_types=np.asarray(section_types, dtype=np.int64),
        soma_points=np.asarray(soma_points, dtype=float).reshape((-1, 3)),
        soma_radii=(np.asarray(soma_diameters, dtype=float) / 2.0),
    )
    tree.validate()
    return tree


def _safe_frame(direction: np.ndarray, parent_rotation: np.ndarray) -> np.ndarray:
    """Parallel-transport-like frame with z following the new edge direction."""
    z = np.asarray(direction, dtype=float)
    n = float(np.linalg.norm(z))
    if n <= 1e-15:
        return parent_rotation.copy()
    z = z / n

    # Preserve the parent's transverse x direction as much as possible.
    ref = parent_rotation[:, 0]
    x = ref - z * float(np.dot(ref, z))
    if np.linalg.norm(x) < 1e-9:
        ref = parent_rotation[:, 1]
        x = ref - z * float(np.dot(ref, z))
    if np.linalg.norm(x) < 1e-9:
        # Deterministic global fallback for the rare parallel degeneracy.
        ref = np.asarray([1.0, 0.0, 0.0])
        if abs(float(np.dot(ref, z))) > 0.9:
            ref = np.asarray([0.0, 1.0, 0.0])
        x = ref - z * float(np.dot(ref, z))

    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    y = y / np.linalg.norm(y)
    # Re-orthogonalize x to suppress accumulated roundoff.
    x = np.cross(y, z)
    x = x / np.linalg.norm(x)
    return np.column_stack([x, y, z])


def build_matrix_scaffold(tree: PointTree) -> MatrixScaffold:
    """Encode a point tree as parent-local SE(3) transforms."""
    tree.validate()
    n = len(tree.positions)
    world_rot = np.zeros((n, 3, 3), dtype=float)
    world_rot[0] = np.eye(3)

    local = np.repeat(np.eye(4, dtype=float)[None, :, :], n, axis=0)

    for i in range(1, n):
        p = int(tree.parents[i])
        Rp = world_rot[p]
        delta = tree.positions[i] - tree.positions[p]
        Rc = _safe_frame(delta, Rp)
        world_rot[i] = Rc

        local[i, :3, :3] = Rp.T @ Rc
        local[i, :3, 3] = Rp.T @ delta

    scaffold = MatrixScaffold(
        root_position=tree.positions[0].copy(),
        parents=tree.parents.copy(),
        local_transforms=local,
        radii=tree.radii.copy(),
        section_ids=tree.section_ids.copy(),
        section_types=tree.section_types.copy(),
    )
    scaffold.validate()
    return scaffold


def reconstruct(scaffold: MatrixScaffold) -> tuple[np.ndarray, np.ndarray]:
    """Reconstruct world points/frames from only root pose + local matrices."""
    scaffold.validate()
    n = len(scaffold.parents)
    world = np.repeat(np.eye(4, dtype=float)[None, :, :], n, axis=0)
    world[0, :3, 3] = scaffold.root_position

    for i in range(1, n):
        p = int(scaffold.parents[i])
        world[i] = world[p] @ scaffold.local_transforms[i]

    return world[:, :3, 3].copy(), world[:, :3, :3].copy()


def child_counts(parents: np.ndarray) -> np.ndarray:
    counts = np.zeros(len(parents), dtype=np.int64)
    for i in range(1, len(parents)):
        counts[int(parents[i])] += 1
    return counts


def descendant_counts(parents: np.ndarray) -> np.ndarray:
    counts = np.zeros(len(parents), dtype=np.int64)
    # Parent-before-child indexing means reverse accumulation is exact.
    for i in range(len(parents) - 1, 0, -1):
        p = int(parents[i])
        counts[p] += counts[i] + 1
    return counts


def descendant_mask(parents: np.ndarray, pivot: int) -> np.ndarray:
    mask = np.zeros(len(parents), dtype=bool)
    mask[pivot] = True
    for i in range(pivot + 1, len(parents)):
        if mask[int(parents[i])]:
            mask[i] = True
    return mask


def choose_twist_pivot(parents: np.ndarray) -> int:
    """Choose a real bifurcation with a sizeable but non-global distal subtree."""
    d = descendant_counts(parents)
    c = child_counts(parents)
    n = len(parents)
    target = max(25.0, 0.15 * n)

    candidates = [
        i for i in range(1, n)
        if c[i] >= 2 and d[i] >= 20 and d[i] <= 0.45 * n
    ]
    if not candidates:
        candidates = [i for i in range(1, n) if d[i] >= 10 and d[i] <= 0.45 * n]
    if not candidates:
        raise ValueError("no suitable internal pivot")

    return min(candidates, key=lambda i: abs(float(d[i]) - target))


def _axis_rotation(axis: str, angle_rad: float) -> np.ndarray:
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    if axis == "x":
        return np.asarray([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)
    if axis == "y":
        return np.asarray([[c, 0, s], [0, 1, 0], [-s, 0, c]], dtype=float)
    if axis == "z":
        return np.asarray([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)
    raise ValueError("axis must be x, y, or z")


def twist_scaffold(
    scaffold: MatrixScaffold,
    pivot: int,
    angle_degrees: float = 20.0,
    axis: str = "y",
) -> MatrixScaffold:
    """Bend/twist one local frame while keeping its attachment point fixed."""
    if pivot <= 0 or pivot >= len(scaffold.parents):
        raise ValueError("pivot must be a non-root node")

    local = scaffold.local_transforms.copy()
    R = _axis_rotation(axis, np.deg2rad(float(angle_degrees)))
    local[pivot, :3, :3] = local[pivot, :3, :3] @ R

    out = MatrixScaffold(
        root_position=scaffold.root_position.copy(),
        parents=scaffold.parents.copy(),
        local_transforms=local,
        radii=scaffold.radii.copy(),
        section_ids=scaffold.section_ids.copy(),
        section_types=scaffold.section_types.copy(),
    )
    out.validate()
    return out


def edge_lengths(positions: np.ndarray, parents: np.ndarray) -> np.ndarray:
    if len(positions) <= 1:
        return np.empty(0, dtype=float)
    out = np.empty(len(positions) - 1, dtype=float)
    for i in range(1, len(positions)):
        out[i - 1] = np.linalg.norm(positions[i] - positions[int(parents[i])])
    return out


def rotation_quality(local_transforms: np.ndarray) -> dict[str, float]:
    errs = []
    dets = []
    I = np.eye(3)
    for T in local_transforms[1:]:
        R = T[:3, :3]
        errs.append(float(np.linalg.norm(R.T @ R - I, ord="fro")))
        dets.append(float(np.linalg.det(R)))
    return {
        "max_orthogonality_error": float(max(errs, default=0.0)),
        "max_abs_det_minus_one": float(max((abs(d - 1.0) for d in dets), default=0.0)),
    }


def export_npz(path: str | Path, scaffold: MatrixScaffold) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        p,
        root_position=scaffold.root_position,
        parents=scaffold.parents,
        local_transforms=scaffold.local_transforms,
        radii=scaffold.radii,
        section_ids=scaffold.section_ids,
        section_types=scaffold.section_types,
    )


def load_npz(path: str | Path) -> MatrixScaffold:
    with np.load(path) as z:
        out = MatrixScaffold(
            root_position=np.asarray(z["root_position"], dtype=float),
            parents=np.asarray(z["parents"], dtype=np.int64),
            local_transforms=np.asarray(z["local_transforms"], dtype=float),
            radii=np.asarray(z["radii"], dtype=float),
            section_ids=np.asarray(z["section_ids"], dtype=np.int64),
            section_types=np.asarray(z["section_types"], dtype=np.int64),
        )
    out.validate()
    return out


def plot_tree(
    path: str | Path,
    positions: np.ndarray,
    parents: np.ndarray,
    *,
    title: str,
    pivot: int | None = None,
) -> None:
    """Fast 3-D line render using one collection rather than one call per edge."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    segs = np.asarray(
        [[positions[int(parents[i])], positions[i]] for i in range(1, len(positions))],
        dtype=float,
    )

    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Line3DCollection(segs, linewidths=0.35))

    mins = np.min(positions, axis=0)
    maxs = np.max(positions, axis=0)
    center = 0.5 * (mins + maxs)
    radius = 0.5 * float(np.max(maxs - mins))
    radius = max(radius, 1.0)
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_zlabel("z (um)")
    ax.set_title(title)
    if pivot is not None:
        ax.scatter(
            [positions[pivot, 0]],
            [positions[pivot, 1]],
            [positions[pivot, 2]],
            s=25,
        )
    fig.tight_layout()
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=180)
    plt.close(fig)
