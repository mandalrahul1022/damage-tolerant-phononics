"""Locked Phase 3 parameters from blueprint Section 3."""

from dataclasses import dataclass


@dataclass(frozen=True)
class LatticeParams:
    a: float = 0.04
    h: float = 0.003
    w_intra: float = 0.0025
    w_inter: float = 0.00375


@dataclass(frozen=True)
class MaterialParams:
    E: float = 3.0e9
    rho: float = 1240.0
    nu: float = 0.36
    eta: float = 0.02


@dataclass(frozen=True)
class SolverParams:
    n_eigs: int = 30
    target_hz: float = 4000.0
    k_path_n: int = 60
