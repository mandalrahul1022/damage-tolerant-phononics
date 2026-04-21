"""Elasticity helpers used across Phase 3 scripts/tests."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt

import numpy as np


@dataclass(frozen=True)
class IsotropicMaterial:
    """Minimal isotropic material model for Phase 3 analytical checks."""

    E: float
    nu: float
    rho: float


def lame_parameters(E: float, nu: float) -> tuple[float, float]:
    """Return Lamé parameters (lambda, mu)."""
    if E <= 0:
        raise ValueError("E must be positive")
    if not (-1.0 < nu < 0.5):
        raise ValueError("nu must be in (-1, 0.5)")

    mu = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return lam, mu


def wave_speeds(E: float, nu: float, rho: float) -> tuple[float, float]:
    """Return longitudinal and transverse speeds (c_L, c_T), m/s."""
    if rho <= 0:
        raise ValueError("rho must be positive")

    lam, mu = lame_parameters(E, nu)
    c_l = sqrt((lam + 2.0 * mu) / rho)
    c_t = sqrt(mu / rho)
    return c_l, c_t


def hz_from_speed_and_wavenumber(speed: float, k_mag: float) -> float:
    """Convert omega = c|k| into Hz."""
    if speed < 0 or k_mag < 0:
        raise ValueError("speed and k_mag must be non-negative")
    return speed * k_mag / (2.0 * pi)


def branch_frequencies_hz(material: IsotropicMaterial, kvec: np.ndarray) -> dict[str, float]:
    """Analytical isotropic elastic branch frequencies for a given wavevector."""
    k = np.asarray(kvec, dtype=float)
    if k.shape != (3,):
        raise ValueError("kvec must be shape (3,)")

    k_mag = float(np.linalg.norm(k))
    c_l, c_t = wave_speeds(material.E, material.nu, material.rho)

    return {
        "transverse_1": hz_from_speed_and_wavenumber(c_t, k_mag),
        "transverse_2": hz_from_speed_and_wavenumber(c_t, k_mag),
        "longitudinal": hz_from_speed_and_wavenumber(c_l, k_mag),
    }
