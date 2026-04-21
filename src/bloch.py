"""Bloch utilities for Phase 3."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def bloch_phase(kvec: np.ndarray, lattice_vec: np.ndarray) -> complex:
    """Return exp(i k·a)."""
    k = np.asarray(kvec, dtype=float)
    a = np.asarray(lattice_vec, dtype=float)
    return complex(np.exp(1j * float(np.dot(k, a))))


@dataclass(frozen=True)
class CubeBlochKPoints:
    gamma: np.ndarray
    x: np.ndarray
    y: np.ndarray


def cube_reference_kpoints(L: float) -> CubeBlochKPoints:
    """Return k=(0,0,0), (pi/L,0,0), (0,pi/L,0)."""
    if L <= 0:
        raise ValueError("L must be positive")

    pi_over_L = np.pi / L
    return CubeBlochKPoints(
        gamma=np.array([0.0, 0.0, 0.0], dtype=float),
        x=np.array([pi_over_L, 0.0, 0.0], dtype=float),
        y=np.array([0.0, pi_over_L, 0.0], dtype=float),
    )


def relative_error(measured: float, expected: float) -> float:
    """Absolute relative error with exact-zero guard."""
    if expected == 0.0:
        return abs(measured)
    return abs(measured - expected) / abs(expected)
