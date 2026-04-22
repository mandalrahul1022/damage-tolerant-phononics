"""Analytical references for Bloch-PBC verification on a homogeneous cube."""

from __future__ import annotations

import numpy as np

from src.params import MaterialParams


def plane_wave_transverse_hz(k_magnitude: float, mat: MaterialParams) -> float:
    """omega_T / (2*pi) where c_T = sqrt(mu / rho), mu = E / (2*(1+nu))."""
    mu = mat.E / (2.0 * (1.0 + mat.nu))
    omega = np.sqrt(mu / mat.rho) * k_magnitude
    return omega / (2.0 * np.pi)


def plane_wave_longitudinal_hz(k_magnitude: float, mat: MaterialParams) -> float:
    """omega_L / (2*pi) where c_L = sqrt((lam + 2*mu) / rho)."""
    mu = mat.E / (2.0 * (1.0 + mat.nu))
    lam = mat.E * mat.nu / ((1.0 + mat.nu) * (1.0 - 2.0 * mat.nu))
    omega = np.sqrt((lam + 2.0 * mu) / mat.rho) * k_magnitude
    return omega / (2.0 * np.pi)
