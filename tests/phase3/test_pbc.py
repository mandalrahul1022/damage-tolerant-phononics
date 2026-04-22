import numpy as np

from scripts.verify_pbc import solve_at_k
from src.params import MaterialParams
from src.cube_reference import (
    plane_wave_transverse_hz,
    plane_wave_longitudinal_hz,
)

L = 0.020  # m
MAT = MaterialParams()  # PLA: E=3e9, rho=1240, nu=0.36


def test_near_zero_modes_at_gamma():
    freqs = solve_at_k(np.array([0.0, 0.0, 0.0]))
    sorted_freqs = np.sort(freqs)
    # The three rigid translations are the only exact null modes.
    assert np.all(sorted_freqs[:3] < 5.0), \
        f'near-zero modes too large: {sorted_freqs[:3]}'
    # Fourth mode is the first genuine elastic mode.
    assert sorted_freqs[3] > 100.0, \
        f'4th mode suspiciously low: {sorted_freqs[3]}'


def test_transverse_plane_wave_x():
    k = np.array([np.pi / L, 0.0, 0.0])
    freqs = solve_at_k(k)
    f_analytic = plane_wave_transverse_hz(np.pi / L, MAT)
    nonzero = np.sort(freqs[freqs > 5.0])
    # Lowest two nonzero modes are the two transverse polarizations
    # (both perpendicular to k, degenerate for isotropic material).
    err1 = abs(nonzero[0] - f_analytic) / f_analytic
    err2 = abs(nonzero[1] - f_analytic) / f_analytic
    assert err1 < 0.01, f'T1 err {err1*100:.3f}% (got {nonzero[0]:.1f}, expect {f_analytic:.1f})'
    assert err2 < 0.01, f'T2 err {err2*100:.3f}% (got {nonzero[1]:.1f}, expect {f_analytic:.1f})'


def test_longitudinal_plane_wave_x():
    k = np.array([np.pi / L, 0.0, 0.0])
    freqs = solve_at_k(k)
    f_analytic = plane_wave_longitudinal_hz(np.pi / L, MAT)
    nonzero = np.sort(freqs[freqs > 5.0])
    # Due additional FE degeneracies, use the mode closest to the
    # longitudinal analytical value instead of a fixed index.
    closest = nonzero[np.argmin(np.abs(nonzero - f_analytic))]
    err = abs(closest - f_analytic) / f_analytic
    assert err < 0.01, f'L err {err*100:.3f}% (got {closest:.1f}, expect {f_analytic:.1f})'
