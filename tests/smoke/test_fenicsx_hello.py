"""
Smoke test: verify FEniCSx (dolfinx) can solve a Poisson problem.

Implements: BLUEPRINT.md Section 5.1 (Phase 2 Step 5.1, Step A4 of the
session plan). Confirms FEniCSx is not just importable but functional.

Verifies:
    Solves -Laplacian(u) = -6 on the unit square with Dirichlet boundary
    condition u(x, y) = 1 + x^2 + 2 y^2 on all of dOmega. The manufactured
    solution is u_exact = 1 + x^2 + 2 y^2 (since -Laplacian of that is -6).
    Uses P2 Lagrange elements on a 32 x 32 mesh; P2 can represent a
    quadratic exactly, so the expected L2 error is near machine precision.

    The test passes if the computed L2 error against the manufactured
    solution is below 1e-3. This is a generous tolerance: we are smoke-
    testing that the full stack (dolfinx + UFL + PETSc LU solve + MPI
    reduction) is wired up correctly, not probing convergence rates.

How to run:
    conda activate kagome-hoti
    python tests/smoke/test_fenicsx_hello.py
or:
    pytest tests/smoke/test_fenicsx_hello.py -v

Expected output on PASS:
    A single line:  "Poisson L2 error: X.XXXe-YY  (tolerance 1e-03)"
    followed by "PASS". Exit code 0. Expected error magnitude on P2 /
    32 x 32 is ~1e-14 (machine precision), well below the 1e-3 threshold.

Pass/fail criteria:
    PASS: solve completes, L2 error < 1e-3, exit 0.
    FAIL: any exception, or L2 error >= 1e-3, or nonzero exit.
"""

from __future__ import annotations

import sys

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, mesh
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI

N_CELLS = 32
DEGREE = 2
TOL = 1e-3


def _poisson_l2_error() -> float:
    msh = mesh.create_unit_square(MPI.COMM_WORLD, N_CELLS, N_CELLS)
    V = fem.functionspace(msh, ("Lagrange", DEGREE))

    def u_exact(x):
        return 1.0 + x[0] ** 2 + 2.0 * x[1] ** 2

    u_D = fem.Function(V)
    u_D.interpolate(u_exact)

    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_connectivity(fdim, tdim)
    boundary_facets = mesh.exterior_facet_indices(msh.topology)
    boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
    bc = fem.dirichletbc(u_D, boundary_dofs)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    f = fem.Constant(msh, default_scalar_type(-6.0))
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = ufl.inner(f, v) * ufl.dx

    problem = LinearProblem(
        a,
        L,
        bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    uh = problem.solve()

    err_form = fem.form(ufl.inner(uh - u_D, uh - u_D) * ufl.dx)
    err_sq = msh.comm.allreduce(fem.assemble_scalar(err_form), op=MPI.SUM)
    return float(np.sqrt(np.real(err_sq)))


def test_poisson_converges() -> None:
    err = _poisson_l2_error()
    assert err < TOL, f"L2 error {err:.3e} exceeds tolerance {TOL:.0e}"


if __name__ == "__main__":
    err = _poisson_l2_error()
    print(f"Poisson L2 error: {err:.3e}  (tolerance {TOL:.0e})")
    if err < TOL:
        print("PASS")
        sys.exit(0)
    print("FAIL")
    sys.exit(1)
