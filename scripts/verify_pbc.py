"""Real Bloch-PBC verification on a homogeneous elastic cube (Track B)."""

from __future__ import annotations

from functools import lru_cache
from typing import Callable

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, mesh
from mpi4py import MPI
from petsc4py import PETSc
from slepc4py import SLEPc

import dolfinx_mpc

from src.params import MaterialParams

L = 0.020
N_HEX = 20
DEGREE = 2
N_EIGS = 15
GEOM_TOL = 1.0e-12


def _lame_parameters(mat: MaterialParams) -> tuple[float, float]:
    mu = mat.E / (2.0 * (1.0 + mat.nu))
    lam = mat.E * mat.nu / ((1.0 + mat.nu) * (1.0 - 2.0 * mat.nu))
    return lam, mu


def _build_cube_space() -> tuple[mesh.Mesh, fem.FunctionSpace]:
    domain = mesh.create_box(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0, 0.0], dtype=np.float64), np.array([L, L, L], dtype=np.float64)],
        [N_HEX, N_HEX, N_HEX],
        cell_type=mesh.CellType.hexahedron,
    )
    V = fem.functionspace(domain, ("Lagrange", DEGREE, (3,)))
    return domain, V


def _make_indicator(mask_fn: Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]) -> Callable[[np.ndarray], np.ndarray]:
    def indicator(x: np.ndarray) -> np.ndarray:
        x0 = np.isclose(x[0], 0.0, atol=GEOM_TOL)
        y0 = np.isclose(x[1], 0.0, atol=GEOM_TOL)
        z0 = np.isclose(x[2], 0.0, atol=GEOM_TOL)
        return mask_fn(x0, y0, z0)

    return indicator


def _relation(dx: float, dy: float, dz: float) -> Callable[[np.ndarray], np.ndarray]:
    def rel(x: np.ndarray) -> np.ndarray:
        out = np.copy(x)
        out[0] += dx
        out[1] += dy
        out[2] += dz
        return out

    return rel


def _build_bloch_mpc(V: fem.FunctionSpace, k_vec: np.ndarray) -> dolfinx_mpc.MultiPointConstraint:
    mpc = dolfinx_mpc.MultiPointConstraint(V)
    phase_x = np.exp(1j * k_vec[0] * L)
    phase_y = np.exp(1j * k_vec[1] * L)
    phase_z = np.exp(1j * k_vec[2] * L)
    bcs: list[fem.DirichletBC] = []

    # Faces (exclude edges/corner so each slave DOF is constrained exactly once).
    mpc.create_periodic_constraint_geometrical(
        V,
        _make_indicator(lambda x0, y0, z0: x0 & (~y0) & (~z0)),
        _relation(L, 0.0, 0.0),
        bcs,
        scale=default_scalar_type(phase_x),
    )
    mpc.create_periodic_constraint_geometrical(
        V,
        _make_indicator(lambda x0, y0, z0: y0 & (~x0) & (~z0)),
        _relation(0.0, L, 0.0),
        bcs,
        scale=default_scalar_type(phase_y),
    )
    mpc.create_periodic_constraint_geometrical(
        V,
        _make_indicator(lambda x0, y0, z0: z0 & (~x0) & (~y0)),
        _relation(0.0, 0.0, L),
        bcs,
        scale=default_scalar_type(phase_z),
    )

    # Edges.
    mpc.create_periodic_constraint_geometrical(
        V,
        _make_indicator(lambda x0, y0, z0: x0 & y0 & (~z0)),
        _relation(L, L, 0.0),
        bcs,
        scale=default_scalar_type(phase_x * phase_y),
    )
    mpc.create_periodic_constraint_geometrical(
        V,
        _make_indicator(lambda x0, y0, z0: x0 & z0 & (~y0)),
        _relation(L, 0.0, L),
        bcs,
        scale=default_scalar_type(phase_x * phase_z),
    )
    mpc.create_periodic_constraint_geometrical(
        V,
        _make_indicator(lambda x0, y0, z0: y0 & z0 & (~x0)),
        _relation(0.0, L, L),
        bcs,
        scale=default_scalar_type(phase_y * phase_z),
    )

    # Corner.
    mpc.create_periodic_constraint_geometrical(
        V,
        _make_indicator(lambda x0, y0, z0: x0 & y0 & z0),
        _relation(L, L, L),
        bcs,
        scale=default_scalar_type(phase_x * phase_y * phase_z),
    )
    mpc.finalize()
    return mpc


@lru_cache(maxsize=1)
def _assemble_once() -> tuple[mesh.Mesh, fem.FunctionSpace, fem.Form, fem.Form]:
    """Build geometry and forms once; only MPC phase changes across k."""
    domain, V = _build_cube_space()
    mat = MaterialParams()
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    lam, mu = _lame_parameters(mat)
    rho = mat.rho
    I = ufl.Identity(3)

    def eps(w):
        return 0.5 * (ufl.grad(w) + ufl.grad(w).T)

    def sigma(w):
        return lam * ufl.tr(eps(w)) * I + 2.0 * mu * eps(w)

    a_K = fem.form(ufl.inner(sigma(u), eps(v)) * ufl.dx)
    a_M = fem.form(rho * ufl.inner(u, v) * ufl.dx)
    return domain, V, a_K, a_M


def _assemble_operators(
    a_K: fem.Form, a_M: fem.Form, mpc: dolfinx_mpc.MultiPointConstraint
) -> tuple[PETSc.Mat, PETSc.Mat]:
    K = dolfinx_mpc.assemble_matrix(a_K, mpc, bcs=[], diagval=1.0)
    # Keep constrained rows out of the generalized mass operator to avoid
    # artificial lambda=1 modes from MPC bookkeeping diagonals.
    M = dolfinx_mpc.assemble_matrix(a_M, mpc, bcs=[], diagval=0.0)
    K.assemble()
    M.assemble()
    return K, M


def _solve_generalized(K: PETSc.Mat, M: PETSc.Mat, n_eigs: int) -> np.ndarray:
    for tol, max_it in ((1.0e-10, 1000), (1.0e-8, 3000)):
        eps = SLEPc.EPS().create(MPI.COMM_WORLD)
        eps.setOperators(K, M)
        eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
        eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
        eps.setDimensions(nev=n_eigs)
        eps.setTolerances(tol=tol, max_it=max_it)
        eps.setTarget(0.0)
        eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)

        st = eps.getST()
        st.setType(SLEPc.ST.Type.SINVERT)
        st.setShift(0.0)
        ksp = st.getKSP()
        ksp.setType(PETSc.KSP.Type.PREONLY)
        pc = ksp.getPC()
        pc.setType(PETSc.PC.Type.LU)
        try:
            pc.setFactorSolverType("mumps")
        except Exception:
            pass

        eps.solve()
        nconv = eps.getConverged()
        if nconv < n_eigs:
            continue

        eigvals_raw = np.array(
            [float(np.real(eps.getEigenvalue(i))) for i in range(n_eigs)]
        )
        lam_scale = max(1.0, abs(eigvals_raw).max())
        eigvals = np.where(
            np.abs(eigvals_raw) < 1.0e-6 * lam_scale,
            0.0,
            eigvals_raw,
        )
        if np.any(eigvals < 0.0):
            raise RuntimeError(
                f"Genuinely negative eigenvalue after clamp: {eigvals}"
            )
        freqs = np.sqrt(eigvals) / (2.0 * np.pi)
        return np.sort(freqs)

    raise RuntimeError(
        f"SLEPc converged fewer than {n_eigs} eigenpairs after retries "
        "(target=0.0)."
    )


@lru_cache(maxsize=16)
def _solve_cached(kx: float, ky: float, kz: float) -> tuple[float, ...]:
    k_vec = np.array([kx, ky, kz], dtype=float)
    _domain, V, a_K, a_M = _assemble_once()
    mpc = _build_bloch_mpc(V, k_vec)
    K, M = _assemble_operators(a_K, a_M, mpc)
    freqs = _solve_generalized(K, M, n_eigs=N_EIGS)
    return tuple(float(x) for x in freqs)


def solve_at_k(k_vec: np.ndarray) -> np.ndarray:
    """Return the 15 lowest Bloch-mode frequencies (Hz) at k-vector `k_vec`."""
    k_arr = np.asarray(k_vec, dtype=float)
    if k_arr.shape != (3,):
        raise ValueError(f"k_vec must have shape (3,), got {k_arr.shape}")
    freqs = _solve_cached(float(k_arr[0]), float(k_arr[1]), float(k_arr[2]))
    return np.asarray(freqs, dtype=float).copy()


if __name__ == "__main__":
    gamma = np.array([0.0, 0.0, 0.0], dtype=float)
    kx = np.array([np.pi / L, 0.0, 0.0], dtype=float)
    print("[k=(0,0,0)]", solve_at_k(gamma))
    print("[k=(pi/L,0,0)]", solve_at_k(kx))
