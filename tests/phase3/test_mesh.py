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


def _assert_periodic_pair(mesh: meshio.Mesh, group_m: str, group_p: str, shift_xy: np.ndarray, tol: float = 1.5e-3) -> None:
    pts_m = _surface_group_points(mesh, group_m)
    pts_p = _surface_group_points(mesh, group_p)

    assert len(pts_m) > 0 and len(pts_p) > 0, f"Empty periodic group in {group_m}/{group_p}"

    shifted = pts_m.copy()
    shifted[:, :2] += shift_xy

    # The + side can be a strict subset depending on clipped half-stubs; enforce
    # that every +side node is represented on the shifted -side.
    tree = cKDTree(shifted)
    dists, _ = tree.query(pts_p, k=1)
    assert float(np.max(dists)) < tol, f"{group_p} is not represented by shifted {group_m}, max_dist={np.max(dists):.3e}"


@pytest.mark.parametrize("ratio", [0.5, 1.0, 1.5])
def test_element_count(ratio: float) -> None:
    mesh = _load_mesh(ratio)
    n_tets = _tet_count(mesh)
    assert 30_000 <= n_tets <= 80_000, f"tets={n_tets} out of [30k, 80k]"


@pytest.mark.parametrize("ratio", [0.5, 1.0, 1.5])
def test_periodic_facet_pairs(ratio: float) -> None:
    mesh = _load_mesh(ratio)
    a1, a2 = lattice_vectors(A)
    _assert_periodic_pair(mesh, "f1m", "f1p", shift_xy=a1)
    _assert_periodic_pair(mesh, "f2m", "f2p", shift_xy=a2)


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
