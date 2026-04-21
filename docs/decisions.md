# Architectural Decision Record

---

## Decision 1: All UFL forms use `ufl.inner`, never `ufl.dot` or `*`

**Date:** 2026-04-21
**Phase:** 2 Step 5.1 (during FEniCSx smoke test)
**Status:** Active, project-wide

### Context
The project environment was deliberately built against PETSc compiled with
complex scalars (`petsc=*=*complex*` in environment.yml) because Phase 4
(viscoelastic damping) requires complex-valued stiffness matrices to
represent E* = E'(1 + i*tan_delta). This is a locked architectural decision
in CLAUDE.md and BLUEPRINT.md Section 3.2.

### Problem encountered
On first run of tests/smoke/test_fenicsx_hello.py, UFL raised
`ArityMismatch: Failure to conjugate test function in complex Form` from
`ufl.algorithms.check_arities`. The offending constructs were `ufl.dot(grad(u),
grad(v))` and `f*v` in the weak form.

### Root cause
When PETSc.ScalarType is complex128, UFL enforces Hermitian inner products
for well-posedness of the variational problem. `ufl.dot` computes a simple
sum product without conjugation; `ufl.inner` applies complex conjugation to
the test function. Real-scalar PETSc builds do not enforce this because
conjugation is a no-op on reals.

### Decision
Every weak form written in this project uses `ufl.inner(...)`, never
`ufl.dot(...)` or bare `*` for vector/tensor contractions with test
functions. This applies to:
  - All Phase 3 Bloch elasticity forms
  - All Phase 4 complex-modulus damping forms
  - All Phase 6 finite-sample eigenvalue forms
  - All future Phase 7 disorder realizations

### Consequences
- Code reviews must flag any `ufl.dot` or bare `*` involving a TestFunction.
- If at some future date we need to temporarily switch to real PETSc (e.g.
  for a performance comparison), this decision must be revisited, not
  silently violated.
- Scalar-valued expressions can still use `*`. This rule applies only where
  a TestFunction or TrialFunction appears in a product.

### Evidence
- tests/smoke/test_fenicsx_hello.py uses ufl.inner throughout after fix,
  now passing with L2 error 2.563e-13 on 32x32 P2 unit square.
- Commit 4bcbbf2.

---
