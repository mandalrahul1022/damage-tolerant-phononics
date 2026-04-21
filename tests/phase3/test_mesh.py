"""Phase 3 Track A gate checks for generated unit-cell meshes."""

from __future__ import annotations

from pathlib import Path

import meshio
import numpy as np
import pytest
from scipy.spatial import cKDTree

from src.mesh_utils import lattice_vectors, rotate_xy, wrap_xy_to_cell

A = 0.04
MESH_DIR = Path("meshes")


def _load_mesh(ratio: float):
    path = MESH_DIR / f"unit_r{ratio}.msh"
    assert path.exists(), f"Missing mesh file: {path}"
    return meshio.read(path)


def _tet_count(mesh: meshio.Mesh) -> int:
    tets = mesh.cells_dict.get("tetra")
    assert tets is not None, "No tetrahedral cells found"
    return int(len(tets))


def _surface_group_points(mesh: meshio.Mesh, group_name: str) -> np.ndarray:
    assert group_name in mesh.field_data, f"Physical group {group_name} not found"
    group_id, _dim = mesh.field_data[group_name]

    tri_blocks = []
    tri_phys = []
    for cells, phys in zip(mesh.cells, mesh.cell_data.get("gmsh:physical", [])):
        if cells.type == "triangle":
            tri_blocks.append(cells.data)
            tri_phys.append(np.asarray(phys, dtype=int))

    assert tri_blocks, "No triangle cell blocks found"

    tri_all = np.concatenate(tri_blocks, axis=0)
    phys_all = np.concatenate(tri_phys, axis=0)
    tri_group = tri_all[phys_all == int(group_id)]
    assert tri_group.size > 0, f"No triangles found for group {group_name}"

    node_ids = np.unique(tri_group.reshape(-1))
    return mesh.points[node_ids, :3]


def _assert_strict_periodic(mesh: meshio.Mesh, src_group: str, dst_group: str, translation_xy: np.ndarray, tol: float = 1e-7) -> None:
    """Bidirectional node-matched periodicity check.

    Asserts that (a) the two physical groups have identical node counts, and
    (b) shifting every src-side node by +translation lands within `tol` of a
    dst-side node, and (c) shifting every dst-side node by -translation lands
    within `tol` of a src-side node. This is the pairing guarantee required
    by dolfinx_mpc.MultiPointConstraint for Bloch PBC assembly.
    """
    src = _surface_group_points(mesh, src_group)
    dst = _surface_group_points(mesh, dst_group)

    assert len(src) == len(dst), (
        f"node count mismatch {src_group}={len(src)} vs {dst_group}={len(dst)}"
    )

    t3 = np.array([translation_xy[0], translation_xy[1], 0.0], dtype=float)

    tree_dst = cKDTree(dst)
    d_fwd, _ = tree_dst.query(src + t3)
    assert float(d_fwd.max()) < tol, (
        f"{src_group}->{dst_group} forward max dist {d_fwd.max():.3e} > tol {tol:.1e}"
    )

    tree_src = cKDTree(src)
    d_rev, _ = tree_src.query(dst - t3)
    assert float(d_rev.max()) < tol, (
        f"{dst_group}->{src_group} reverse max dist {d_rev.max():.3e} > tol {tol:.1e}"
    )


@pytest.mark.parametrize("ratio", [0.5, 1.0, 1.5])
def test_element_count(ratio: float) -> None:
    mesh = _load_mesh(ratio)
    n_tets = _tet_count(mesh)
    assert 30_000 <= n_tets <= 80_000, f"tets={n_tets} out of [30k, 80k]"


@pytest.mark.parametrize("ratio", [0.5, 1.0, 1.5])
def test_exact_node_periodicity_a1(ratio: float) -> None:
    mesh = _load_mesh(ratio)
    a1, _ = lattice_vectors(A)
    _assert_strict_periodic(mesh, "f1m", "f1p", translation_xy=a1)


@pytest.mark.parametrize("ratio", [0.5, 1.0, 1.5])
def test_exact_node_periodicity_a2(ratio: float) -> None:
    mesh = _load_mesh(ratio)
    _, a2 = lattice_vectors(A)
    _assert_strict_periodic(mesh, "f2m", "f2p", translation_xy=a2)


def test_c3_symmetry_at_ratio_1() -> None:
    mesh = _load_mesh(1.0)
    points = mesh.points[:, :3]

    a1, a2 = lattice_vectors(A)
    centroid = (a1 + a2) / 6.0

    # Check symmetry on the core motif region, where the breathing kagome
    # structure is expected to satisfy C3 at ratio=1.
    radius = np.linalg.norm(points[:, :2] - centroid, axis=1)
    core = points[radius < 0.012]

    wrapped = core.copy()
    wrapped[:, :2] = wrap_xy_to_cell(wrapped[:, :2], a=A)

    rotated = core.copy()
    rotated[:, :2] = rotate_xy(core[:, :2], angle_rad=2.0 * np.pi / 3.0, center=centroid)
    rotated[:, :2] = wrap_xy_to_cell(rotated[:, :2], a=A)

    tree = cKDTree(wrapped)
    dists, _ = tree.query(rotated, k=1)
    assert float(np.percentile(dists, 95.0)) < 2.0e-3, (
        f"C3 wrapped-rotation mismatch, p95_dist={np.percentile(dists, 95.0):.3e}"
    )
