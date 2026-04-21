# Kagome HOTI Simulation

A 25-week computational physics study of a printable audio-frequency
**breathing kagome higher-order topological insulator plate**, modeled as a
fully 3D, viscoelastically-damped, orthotropic elastic system. The project
produces five quantitative predictions suitable for submission as a
first-authored paper to *Physical Review Applied*.

Simulation-only. No fabrication.

## Advisor

Dr. Yang, Caldwell.

## Setup

```bash
conda env create -f environment.yml
conda activate kagome-hoti
```

Verify FEniCSx is importable:

```bash
python -c "import dolfinx; print(dolfinx.__version__)"
```

## Layout

- `CLAUDE.md` — project-wide context loaded into every Claude Code session
- `BLUEPRINT.md` — authoritative 27-page project plan
- `STATUS.md` — current phase, atomic step, last verified checkpoint
- `src/` — library code
- `scripts/` — runnable entry points (Bloch solves, sweeps, figures)
- `meshes/` — Gmsh-generated `.msh`/`.xdmf` (not tracked)
- `results/` — numerical outputs (not tracked)
- `figures/` — publication-ready figures (not tracked)
- `notebooks/` — exploratory notebooks + `log.md` lab notebook
- `tests/` — correctness tests; `tests/smoke/` are fast phase-gate checks
- `docs/theory.md` — derivation traces
- `docs/decisions.md` — architectural decision record
- `.claude/commands/` — project slash commands (`/start-session`, `/verify-phase`, `/end-session`)

## Citation

> Mandal, R. R. *Higher-order topological corner states in a viscoelastically
> damped breathing kagome plate: a 3D computational study.* In preparation,
> 2026. (Target venue: Physical Review Applied.)
