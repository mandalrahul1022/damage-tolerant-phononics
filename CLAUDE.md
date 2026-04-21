# Kagome HOTI Simulation - Project Context for Claude

## Who you are helping
Rahul Mandal, undergraduate researcher at Caldwell, starting this project
from ground zero at ~5 hours per week. Faculty advisor: Dr. Yang. Target venue:
Physical Review Applied. Timeline: 25 weeks to submission.

## What we are building
A fully 3D, viscoelastically-damped, orthotropic computational framework for a
printable audio-frequency breathing kagome higher-order topological insulator
plate, producing five specific quantitative predictions publishable as a
standalone computational physics paper.

## Core physics you must respect
- The lattice is BREATHING KAGOME, not honeycomb valley Hall. Do not propose
  valley Chern number anywhere. The valid topological invariant is the Z3 Berry
  phase computed via discretized Wilson loop.
- The protection is crystalline C3 symmetry, not valley conservation.
- Corner states appear at 60-degree ACUTE terminations only. 120-degree obtuse
  corners give zero corner states. Finite samples are equilateral triangles.
- The project is SIMULATION ONLY. No physical fabrication work. Do not propose
  experimental protocols, hardware, or printing steps.
- Damping matters. Kelvin-Voigt is a reference baseline; Zener SLS is the real
  model. Expect non-monotonic gap-vs-damping under Zener (branch overtaking).
- The simulation is 3D VECTOR elasticity. 2D plane-stress misses out-of-plane
  extensional and shear modes that can fill the flexural gap experimentally.

## Locked parameters (do not change without flagging)
- Lattice constant a = 4 cm
- Plate thickness h = 3 mm
- Nominal breathing ratio = 1.5 (sweep 0.5 to 2.0)
- Intra-cell beam width = 2.5 mm, inter-cell = 3.75 mm
- Material: PLA, E' = 3.0 GPa, rho = 1240 kg/m^3, nu = 0.36
- Nominal loss tangent = 0.02 (sweep 0 to 0.08)
- Finite sample: equilateral triangle, 12 unit cells per side
- Target gap center: 3.5-4.5 kHz

## The toolchain
- Python 3.11 in a conda env named `kagome-hoti`
- FEniCSx (dolfinx) for finite element Bloch solves
- Gmsh Python API for parametric meshes
- SLEPc for eigenvalue problems (complex-valued for damping phase)
- NumPy/SciPy/matplotlib for everything else
- PETSc under the hood, CPU-bound (GPU does not help FEniCSx Bloch solves)

## Current phase
See STATUS.md for current phase and atomic step. Do not start work until you
have read STATUS.md.

## Rules for writing code in this repo
- Every atomic step in BLUEPRINT.md must be verifiable. If you write code, also
  write the smoke test that verifies it.
- Cross-check every FEniCSx result against the NumPy tight-binding reference
  in tests/tb_reference.py before trusting it.
- Commit after every verified atomic step. Commit messages reference the
  blueprint section number (e.g. "Phase 3 Step 6.2.4: assemble K and M at Gamma").
- Do not move to the next phase until ALL verification criteria for the
  current phase are checked off. Phase gates are hard gates.
- If a verification criterion fails, first consult the Failure Modes subsection
  for that phase in BLUEPRINT.md. 80% of issues are already listed there.

## What to escalate to me (Raj), not resolve yourself
- Any change to locked parameters in the list above
- Any decision to skip or defer a verification criterion
- Any physics interpretation question where published literature is ambiguous
- When a phase gate fails and the blueprint's failure modes do not match the
  observed problem

## Things you should do automatically
- Read STATUS.md before responding
- Suggest which blueprint section is relevant to the current question
- Run smoke tests after code edits before reporting success
- Update STATUS.md at session end if I ask you to

## Reference documents
- BLUEPRINT.md: full project plan, 27 pages, authoritative
- STATUS.md: current state, updated each session
- docs/theory.md: derivation traces I have completed
- docs/decisions.md: architectural decisions made so far
- notebooks/log.md: my running lab notebook