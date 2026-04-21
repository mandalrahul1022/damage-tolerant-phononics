# Project Status

**Project:** Breathing Kagome HOTI - Computational Blueprint
**PI:** Raj Rahul Mandal  |  **Advisor:** Dr. Yang
**Target:** Physical Review Applied
**Blueprint:** `BLUEPRINT.md`

## Current phase

Phase 2 (Environment Setup and Tight-Binding Sanity Check) - in progress.
Blueprint Section 5, Week 1.

## Phase 2 gate status

| Criterion | Status | Notes |
|---|---|---|
| Conda env activates (local + workstation) | PASS (local); PARTIAL (workstation untested this session) | `environment.yml` committed and local checks pass. Campus workstation verification pending next on-site session. |
| FEniCSx hello world | PASS | `tests/smoke/test_fenicsx_hello.py` |
| Gmsh 2D mesh smoke test | PASS | `tests/smoke/test_gmsh_hello.py` |
| Dirac cone closure at r=1 | PASS | `tests/tb_reference.py`, Criterion 1 |
| Symmetric gap at r=0.5 and r=1.5 | PASS | Canonical Bloch Hamiltonian applied; bulk gaps now match to machine precision (`1.500000`, `1.500000`). |
| Berry phase = 0 at r=0.5 (trivial) | PASS | Wilson loop (`N=50`) |
| Berry phase = 2*pi/3 at r=1.5 (topological) | PASS | Wilson loop (`N=50`) |
| Repo initialized with scaffold | PASS | `src/`, `scripts/`, `meshes/`, `results/`, `figures/`, `notebooks/`, `tests/`, `docs/` |

## Next actions

1. Verify the same Phase 2 checks on the campus GPU workstation.
2. Commit Phase 2 TB fixes and regenerated Phase 2 figures.
3. Start Phase 3 Section 6.1: parametric Gmsh 3D unit-cell mesh script.

## Known issues / untracked work

- `tests/tb_reference.py` and `figures/phase2/*` are currently untracked.
- `node_modules/` and `package*.json` were unrelated to this project and are now ignored.

## Hours log (Phase 2)

Week 1: ~5 hours (environment + smoke tests + TB reference corrections).
