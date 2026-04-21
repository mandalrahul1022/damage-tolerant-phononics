"""
Smoke test: verify the kagome-hoti conda environment is functional.

Implements: BLUEPRINT.md Section 5.1 (Phase 2 Step 5.1, Step A3 of the
session plan) -- post-env-create verification.

Verifies: every load-bearing scientific Python package imports cleanly
and reports a non-empty __version__. This is an import-and-introspect
test only; it does not exercise numerical capability. Numerical capability
is checked by test_fenicsx_hello.py (FEniCSx Poisson) and
test_gmsh_hello.py (Gmsh mesh).

How to run:
    conda activate kagome-hoti
    python tests/smoke/test_env.py
or:
    pytest tests/smoke/test_env.py -v

Expected output on PASS:
    One line per package (dolfinx, gmsh, mpi4py, petsc4py, slepc4py,
    numpy, scipy, matplotlib) printing the package name and its version,
    followed by "ALL IMPORTS OK". Exit code 0.

Pass/fail criteria:
    PASS: every import succeeds, every __version__ is non-empty, exit 0.
    FAIL: any ImportError, missing __version__ attribute, or nonzero exit.
"""

from __future__ import annotations

import sys


def test_imports() -> None:
    import dolfinx
    import gmsh
    import matplotlib
    import mpi4py
    import numpy
    import petsc4py
    import scipy
    import slepc4py

    versions = {
        "dolfinx": dolfinx.__version__,
        "gmsh": gmsh.__version__,
        "mpi4py": mpi4py.__version__,
        "petsc4py": petsc4py.__version__,
        "slepc4py": slepc4py.__version__,
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }

    for name, version in versions.items():
        print(f"  {name:12s} {version}")
        assert version, f"{name} has empty __version__"


if __name__ == "__main__":
    print("kagome-hoti environment smoke test")
    print(f"  {'python':12s} {sys.version.split()[0]}")
    try:
        test_imports()
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        sys.exit(1)
    print("ALL IMPORTS OK")
    sys.exit(0)
