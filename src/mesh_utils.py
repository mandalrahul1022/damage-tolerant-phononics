"""Geometry helpers for Phase 3 unit-cell meshing/tests."""

from __future__ import annotations

import numpy as np


def lattice_vectors(a: float) -> tuple[np.ndarray, np.ndarray]:
    a1 = np.array([a, 0.0], dtype=float)
    a2 = np.array([0.5 * a, np.sqrt(3.0) * 0.5 * a], dtype=float)
    return a1, a2


def lattice_matrix(a: float) -> np.ndarray:
    a1, a2 = lattice_vectors(a)
    return np.column_stack((a1, a2))


def barycentric_coords(xy: np.ndarray, a: float) -> np.ndarray:
    """Map XY coordinates to coefficients (s, t) in xy = s*a1 + t*a2."""
    mat_inv = np.linalg.inv(lattice_matrix(a))
    return np.asarray(xy, dtype=float) @ mat_inv.T


def xy_from_bary(st: np.ndarray, a: float) -> np.ndarray:
    """Map lattice coefficients (s, t) back to Cartesian XY."""
    return np.asarray(st, dtype=float) @ lattice_matrix(a).T


def wrap_xy_to_cell(xy: np.ndarray, a: float) -> np.ndarray:
    """Wrap XY coordinates into the primitive cell [0,1) x [0,1) in (s,t)."""
    st = barycentric_coords(xy, a)
    st_wrapped = np.mod(st, 1.0)
    return xy_from_bary(st_wrapped, a)


def rotate_xy(xy: np.ndarray, angle_rad: float, center: np.ndarray) -> np.ndarray:
    """Rotate XY points around a center."""
    pts = np.asarray(xy, dtype=float) - np.asarray(center, dtype=float)
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    rot = np.array([[c, -s], [s, c]], dtype=float)
    return pts @ rot.T + center
