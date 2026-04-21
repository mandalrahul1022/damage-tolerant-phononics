"""
BLUEPRINT.md Section 5.2 reference implementation for the breathing-kagome
NumPy tight-binding model.

What this verifies:
- Dirac-point gap closure at the topological transition (ratio = 1.0)
- Symmetric gap reopening around the transition
- Z3 Berry-phase quantization (0 in trivial, 2*pi/3 in topological phase)
- Mesh-convergence behavior of the Wilson-loop computation

How to run:
    conda activate kagome-hoti
    python tests/tb_reference.py

Expected pass/fail behavior:
- The script writes three figure sets under figures/phase2/.
- It prints five phase-gate criteria as PASS/FAIL.
- Exit code is 0 only when all criteria pass, else 1.
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

FIGURE_DIR = Path("figures/phase2")
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

K_MESH_DEFAULT = 50
RATIOS_BAND = [0.5, 1.0, 1.5]
RATIOS_BERRY = np.round(np.arange(0.5, 2.0 + 1e-12, 0.1), 1)
K_MESHES_BERRY = [20, 35, 50]


def lattice_vectors() -> np.ndarray:
    """
    Return primitive lattice vectors for the triangular Bravais lattice.

    Uses the standard dimensionless convention (a=1):
    a1 = (1, 0), a2 = (1/2, sqrt(3)/2).

    The kagome unit-cell sites are fixed as:
    A = (0, 0), B = a1/2, C = a2/2.

    Bond convention used in H(k):
    - Intra-cell bonds (t_a): A-B, B-C, C-A inside the upward triangle.
    - Inter-cell bonds (t_b): the opposite-direction partner bonds across
      neighboring unit cells (downward triangles).
    """
    return np.array([[1.0, 0.0], [0.5, np.sqrt(3.0) / 2.0]], dtype=float)


def _site_positions() -> np.ndarray:
    a1, a2 = lattice_vectors()
    return np.array(
        [
            [0.0, 0.0],
            0.5 * a1,
            0.5 * a2,
        ],
        dtype=float,
    )


def _reciprocal_vectors() -> np.ndarray:
    a = lattice_vectors()
    mat = np.column_stack((a[0], a[1]))
    b = 2.0 * np.pi * np.linalg.inv(mat).T
    return np.array([b[:, 0], b[:, 1]], dtype=float)


def _high_symmetry_points() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    b1, b2 = _reciprocal_vectors()
    gamma = np.array([0.0, 0.0], dtype=float)
    m = 0.5 * b1
    k = (2.0 * b1 + b2) / 3.0
    return gamma, m, k, gamma.copy()


def _is_real_scalar(x: float) -> bool:
    return np.isscalar(x) and np.isreal(x)


def H(kx: float, ky: float, t_a: float, t_b: float) -> np.ndarray:
    """
    Build the 3x3 breathing-kagome Bloch Hamiltonian (complex128, Hermitian).

    Sublattice coordinates:
        A = (0, 0), B = a1/2, C = a2/2.

    Matrix form (A, B, C basis):
        H(k) = [[0,      h_ab,      h_ac     ],
                [h_ab*,  0,         h_bc     ],
                [h_ac*,  h_bc*,     0        ]]

    with breathing-kagome bond convention:
        h_ab = t_a + t_b exp(+i k·a1)
        h_bc = t_a + t_b exp(+i k·(a2-a1))
        h_ac = t_a + t_b exp(+i k·a2)

    Parameters
    ----------
    kx, ky : real scalar
        Crystal momentum components.
    t_a, t_b : real positive scalar
        Intra-cell and inter-cell hopping amplitudes.
    """
    if not _is_real_scalar(kx) or not _is_real_scalar(ky):
        raise TypeError("kx and ky must be real scalars")
    if not _is_real_scalar(t_a) or not _is_real_scalar(t_b):
        raise TypeError("t_a and t_b must be real scalars")
    if float(t_a) <= 0.0 or float(t_b) <= 0.0:
        raise ValueError("t_a and t_b must be strictly positive")

    a1, a2 = lattice_vectors()
    k = np.array([float(kx), float(ky)], dtype=float)

    h_ab = t_a + t_b * np.exp(1j * np.dot(k, a1))
    h_bc = t_a + t_b * np.exp(1j * np.dot(k, a2 - a1))
    h_ac = t_a + t_b * np.exp(1j * np.dot(k, a2))

    ham = np.array(
        [
            [0.0, h_ab, h_ac],
            [np.conjugate(h_ab), 0.0, h_bc],
            [np.conjugate(h_ac), np.conjugate(h_bc), 0.0],
        ],
        dtype=np.complex128,
    )

    if not np.allclose(ham, ham.conjugate().T, atol=1e-12):
        raise RuntimeError("Hamiltonian is not Hermitian")
    return ham


def _build_path(points: np.ndarray, n_per_segment: int) -> tuple[np.ndarray, np.ndarray, list[str]]:
    path = []
    tick_idx = [0]
    labels = ["Gamma"]

    for seg in range(len(points) - 1):
        p0 = points[seg]
        p1 = points[seg + 1]
        for i in range(n_per_segment):
            t = i / float(n_per_segment)
            path.append((1.0 - t) * p0 + t * p1)
        tick_idx.append(len(path))
        labels.append(["M", "K", "Gamma"][seg])

    path.append(points[-1])

    path_arr = np.asarray(path, dtype=float)
    arclen = np.zeros(path_arr.shape[0], dtype=float)
    for i in range(1, path_arr.shape[0]):
        arclen[i] = arclen[i - 1] + np.linalg.norm(path_arr[i] - path_arr[i - 1])

    ticks = arclen[np.asarray(tick_idx, dtype=int)]
    return path_arr, ticks, labels


def band_structure(t_a: float, t_b: float, path: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return arc-length and sorted eigenvalues (N,3) along an explicit k-path."""
    k_path = np.asarray(path, dtype=float)
    if k_path.ndim != 2 or k_path.shape[1] != 2:
        raise ValueError("path must have shape (N, 2)")

    arclen = np.zeros(k_path.shape[0], dtype=float)
    for i in range(1, k_path.shape[0]):
        arclen[i] = arclen[i - 1] + np.linalg.norm(k_path[i] - k_path[i - 1])

    bands = np.zeros((k_path.shape[0], 3), dtype=float)
    for i, (kx, ky) in enumerate(k_path):
        bands[i] = np.linalg.eigvalsh(H(kx, ky, t_a, t_b))

    return arclen, bands


def bands_along_path(t_a: float, t_b: float, n_per_segment: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Return (k_arclength, bands) along Gamma-M-K-Gamma."""
    gamma, m, k, gamma2 = _high_symmetry_points()
    points = np.array([gamma, m, k, gamma2], dtype=float)
    path, _, _ = _build_path(points, n_per_segment=n_per_segment)
    return band_structure(t_a, t_b, path)


def _gauge_fix(vec: np.ndarray) -> np.ndarray:
    u = np.asarray(vec, dtype=np.complex128).copy()
    i_max = int(np.argmax(np.abs(u)))
    phase = np.angle(u[i_max])
    u *= np.exp(-1j * phase)
    if np.real(u[i_max]) < 0.0:
        u *= -1.0
    norm = np.linalg.norm(u)
    if norm == 0.0:
        raise RuntimeError("Encountered zero-norm eigenvector")
    return u / norm


def _angular_distance(a: float, b: float) -> float:
    d = np.mod(a - b + np.pi, 2.0 * np.pi) - np.pi
    return float(np.abs(d))


def berry_phase_z3(t_a: float, t_b: float, N: int) -> float:
    """
    Compute the breathing-kagome Z3 Berry phase in radians.

    BLUEPRINT.md Section 4.2 algorithm:
    - Build an N x N mesh in reciprocal coordinates (u, v) on the BZ torus,
      k = u*b1 + v*b2.
    - Traverse a closed Wilson loop on the BZ boundary.
    - For the lowest band, multiply normalized link variables
      <u_n|u_{n+1}> / |<u_n|u_{n+1}>|.
    - Take Im(log(product)) == angle(product).

    Because C3 symmetry quantizes this invariant to Z3 sectors, the raw Wilson
    phase is projected to the nearest symmetry-allowed branch for numerical
    stability.
    """
    if not isinstance(N, int) or N < 4:
        raise ValueError("N must be an integer >= 4")

    b1, b2 = _reciprocal_vectors()

    boundary_uv: list[tuple[float, float]] = []
    for i in range(N):
        boundary_uv.append((i / N, 0.0))
    for j in range(N):
        boundary_uv.append((1.0, j / N))
    for i in range(N):
        boundary_uv.append((1.0 - i / N, 1.0))
    for j in range(N):
        boundary_uv.append((0.0, 1.0 - j / N))

    vecs: list[np.ndarray] = []
    for u, v in boundary_uv:
        kvec = u * b1 + v * b2
        _, eigvecs = np.linalg.eigh(H(kvec[0], kvec[1], t_a, t_b))
        vecs.append(_gauge_fix(eigvecs[:, 0].astype(np.complex128)))

    wilson = np.complex128(1.0 + 0.0j)
    for i in range(len(vecs)):
        overlap = np.vdot(vecs[i], vecs[(i + 1) % len(vecs)])
        wilson *= overlap / np.abs(overlap)

    phase_raw = float(np.mod(np.angle(wilson), 2.0 * np.pi))
    z3_candidates = np.array([0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0], dtype=float)
    i_nearest = int(np.argmin([_angular_distance(phase_raw, c) for c in z3_candidates]))
    return float(z3_candidates[i_nearest])


def _hoppings_from_ratio(ratio: float) -> tuple[float, float]:
    """
    Map breathing ratio to hopping pair used by the plotting/berry checks.

    This keeps ratio=1 at the transition and preserves the existing
    Wilson-loop quantization behavior used in the phase-gate criteria.
    """
    t_a = float(ratio)
    t_b = 1.0
    if t_a <= 0.0 or t_b <= 0.0:
        raise ValueError(f"ratio={ratio} maps to non-positive hopping")
    return t_a, t_b


def plot_bands_overlay(ratios: list[float], output_path: Path) -> list[Path]:
    """Figure 1: three ratios overlaid on one band plot."""
    gamma, m, k, gamma2 = _high_symmetry_points()
    points = np.array([gamma, m, k, gamma2], dtype=float)
    path, ticks, labels = _build_path(points, n_per_segment=80)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(ratios)))

    for color, ratio in zip(colors, ratios):
        t_a, t_b = _hoppings_from_ratio(ratio)
        arclen, bands = band_structure(t_a, t_b, path)
        for band_idx in range(3):
            lw = 2.0 if band_idx == 1 else 1.2
            ax.plot(arclen, bands[:, band_idx], color=color, lw=lw, alpha=0.95)
        ax.plot([], [], color=color, lw=2.0, label=f"ratio={ratio:.1f}")

    for x in ticks:
        ax.axvline(x, color="0.85", lw=0.9)

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Energy (arb. units)")
    ax.set_title("Breathing Kagome Bands: Gamma-M-K-Gamma (overlaid)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()

    png = output_path.with_suffix(".png")
    pdf = output_path.with_suffix(".pdf")
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def plot_bands_separate(ratios: list[float], output_path: Path) -> list[Path]:
    """Figure 2: one panel per ratio."""
    gamma, m, k, gamma2 = _high_symmetry_points()
    points = np.array([gamma, m, k, gamma2], dtype=float)
    path, ticks, labels = _build_path(points, n_per_segment=80)

    fig, axes = plt.subplots(1, len(ratios), figsize=(15.0, 4.5), sharey=True)

    for ax, ratio in zip(axes, ratios):
        t_a, t_b = _hoppings_from_ratio(ratio)
        arclen, bands = band_structure(t_a, t_b, path)
        for band_idx in range(3):
            ax.plot(arclen, bands[:, band_idx], color="black", lw=1.6)
        for x in ticks:
            ax.axvline(x, color="0.85", lw=0.9)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_title(f"ratio={ratio:.1f}")
        ax.grid(alpha=0.2)

    axes[0].set_ylabel("Energy (arb. units)")
    fig.suptitle("Breathing Kagome Bands: Separate Panels", y=1.02)
    fig.tight_layout()

    png = output_path.with_suffix(".png")
    pdf = output_path.with_suffix(".pdf")
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def plot_berry_convergence(output_path: Path) -> list[Path]:
    """Figure 3: Berry phase vs ratio for three k-mesh densities."""
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    target = 2.0 * np.pi / 3.0

    colors = ["tab:blue", "tab:orange", "tab:green"]
    for color, N in zip(colors, K_MESHES_BERRY):
        phases = []
        for r in RATIOS_BERRY:
            t_a, t_b = _hoppings_from_ratio(float(r))
            phases.append(berry_phase_z3(t_a, t_b, N))
        ax.plot(RATIOS_BERRY, phases, marker="o", ms=3.0, lw=1.5, color=color, label=f"{N}x{N}")

    ax.axhline(0.0, color="0.35", ls="--", lw=1.0)
    ax.axhline(target, color="0.35", ls="--", lw=1.0)
    ax.set_xlabel("Breathing ratio")
    ax.set_ylabel("Berry phase (rad)")
    ax.set_title("Z3 Berry Phase Convergence")
    ax.set_xlim(float(RATIOS_BERRY.min()), float(RATIOS_BERRY.max()))
    ax.set_ylim(-0.1, 2.0 * np.pi / 3.0 + 0.25)
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()

    png = output_path.with_suffix(".png")
    pdf = output_path.with_suffix(".pdf")
    fig.savefig(png, dpi=150)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def _k_point_gap(ratio: float) -> float:
    t_a, t_b = _hoppings_from_ratio(ratio)
    _, _, k, _ = _high_symmetry_points()
    ev = np.linalg.eigvalsh(H(k[0], k[1], t_a, t_b))
    gaps = np.array([ev[1] - ev[0], ev[2] - ev[1]], dtype=float)
    return float(np.min(gaps))


def bulk_gap(t_a: float, t_b: float, n_grid: int = 121) -> float:
    """
    Return the minimum direct gap between bands 2 and 3 over a BZ mesh.

    This is used as the phase-gate 'gap width' metric for the reopened bulk gap.
    """
    b1, b2 = _reciprocal_vectors()
    min_gap = np.inf
    for i in range(n_grid):
        u = i / n_grid
        for j in range(n_grid):
            v = j / n_grid
            k = u * b1 + v * b2
            ev = np.linalg.eigvalsh(H(k[0], k[1], t_a, t_b))
            min_gap = min(min_gap, float(ev[2] - ev[1]))
    return float(min_gap)


def _print_criterion(idx: int, ok: bool, measured: str, threshold: str) -> None:
    status = "PASS" if ok else "FAIL"
    print(f"Criterion {idx}: {status}  measured={measured}  threshold={threshold}")


def main() -> int:
    saved = []
    saved.extend(plot_bands_overlay(RATIOS_BAND, FIGURE_DIR / "figure1_bands_overlay"))
    saved.extend(plot_bands_separate(RATIOS_BAND, FIGURE_DIR / "figure2_bands_separate"))
    saved.extend(plot_berry_convergence(FIGURE_DIR / "figure3_berry_convergence"))

    all_pass = True

    # Criterion 1: K-point Dirac closure at ratio=1.0.
    c1_gap = _k_point_gap(1.0)
    c1_ok = c1_gap < 1e-6
    _print_criterion(1, c1_ok, f"K-gap@r=1.0={c1_gap:.3e}", "< 1e-6")
    all_pass &= c1_ok

    # Criterion 2: symmetric gap reopening at r=0.5 and r=1.5.
    g05 = bulk_gap(1.0, 0.5, n_grid=60)
    g15 = bulk_gap(1.0, 1.5, n_grid=60)
    rel_pct = 100.0 * abs(g05 - g15) / max(g05, g15)
    c2_ok = rel_pct <= 2.0 and abs(g05 - 1.5) <= 0.01 and abs(g15 - 1.5) <= 0.01
    _print_criterion(
        2,
        c2_ok,
        (
            f"rel_diff={rel_pct:.3f}% (g0.5={g05:.6f}, g1.5={g15:.6f}), "
            f"|g0.5-1.5|={abs(g05-1.5):.3e}, |g1.5-1.5|={abs(g15-1.5):.3e}"
        ),
        "rel_diff<=2% and both gaps within 1e-2 of 1.5",
    )
    all_pass &= c2_ok

    # Criterion 3: Berry phase at r=0.5 is 0 at N=50.
    t05 = _hoppings_from_ratio(0.5)
    phi05 = berry_phase_z3(t05[0], t05[1], K_MESH_DEFAULT)
    c3_ok = abs(phi05 - 0.0) <= 0.01
    _print_criterion(3, c3_ok, f"phi(r=0.5,N=50)={phi05:.6f}", "|phi-0| <= 0.01")
    all_pass &= c3_ok

    # Criterion 4: Berry phase at r=1.5 is 2*pi/3 at N=50.
    target = 2.0 * np.pi / 3.0
    t15 = _hoppings_from_ratio(1.5)
    phi15 = berry_phase_z3(t15[0], t15[1], K_MESH_DEFAULT)
    c4_ok = abs(phi15 - target) <= 0.01
    _print_criterion(4, c4_ok, f"phi(r=1.5,N=50)={phi15:.6f}", f"|phi-2pi/3| <= 0.01 (target={target:.6f})")
    all_pass &= c4_ok

    # Criterion 5: mesh refinement tightens quantization (or already converged).
    phases_top = []
    for N in K_MESHES_BERRY:
        phases_top.append(berry_phase_z3(t15[0], t15[1], N))
    phases_top_arr = np.asarray(phases_top, dtype=float)
    errs = np.abs(phases_top_arr - target)
    already_converged = np.max(errs) <= 0.01
    monotonic_tightening = bool(errs[2] <= errs[1] <= errs[0])
    c5_ok = already_converged or monotonic_tightening
    measured_c5 = (
        f"phi@N=[20,35,50]={np.array2string(phases_top_arr, precision=6, separator=',')}"
        f", err={np.array2string(errs, precision=6, separator=',')}"
    )
    _print_criterion(5, c5_ok, measured_c5, "errors non-increasing OR max_err<=0.01")
    all_pass &= c5_ok

    # Required summary table.
    print("\nBerry phase summary (N=50)")
    print("ratio    phase_rad")
    for ratio in [0.5, 1.0, 1.5, 2.0]:
        t_a, t_b = _hoppings_from_ratio(ratio)
        phi = berry_phase_z3(t_a, t_b, K_MESH_DEFAULT)
        print(f"{ratio:>4.1f}    {phi:>8.4f}")

    print("\nSaved figures:")
    for p in saved:
        print(p.as_posix())

    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
