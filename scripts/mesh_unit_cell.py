#!/usr/bin/env python3
"""Phase 3 Track A: parametric 3D unit-cell mesh generator."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import gmsh
import meshio
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.mesh_utils import barycentric_coords, lattice_vectors
from src.params import LatticeParams


def _add_cell_prism(a: float, h: float, lc: float) -> int:
    a1, a2 = lattice_vectors(a)
    p0 = gmsh.model.occ.addPoint(0.0, 0.0, 0.0, lc)
    p1 = gmsh.model.occ.addPoint(a1[0], a1[1], 0.0, lc)
    p2 = gmsh.model.occ.addPoint(a1[0] + a2[0], a1[1] + a2[1], 0.0, lc)
    p3 = gmsh.model.occ.addPoint(a2[0], a2[1], 0.0, lc)

    l0 = gmsh.model.occ.addLine(p0, p1)
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p0)
    loop = gmsh.model.occ.addCurveLoop([l0, l1, l2, l3])
    surf = gmsh.model.occ.addPlaneSurface([loop])
    ext = gmsh.model.occ.extrude([(2, surf)], 0.0, 0.0, h)
    vols = [tag for dim, tag in ext if dim == 3]
    return vols[0]


def _beam_corners(start: np.ndarray, end: np.ndarray, width: float) -> np.ndarray:
    d = end - start
    L = float(np.linalg.norm(d))
    if L <= 1e-12:
        raise ValueError("Beam endpoints must be distinct")

    t = d / L
    n = np.array([-t[1], t[0]], dtype=float)
    o = 0.5 * width * n

    return np.array([start + o, end + o, end - o, start - o], dtype=float)


def _add_beam_volume(start: np.ndarray, end: np.ndarray, width: float, h: float, lc: float) -> tuple[int, list[int]]:
    corners = _beam_corners(start=start, end=end, width=width)
    point_tags = [gmsh.model.occ.addPoint(float(p[0]), float(p[1]), 0.0, lc) for p in corners]

    line_tags = []
    for i in range(4):
        line_tags.append(gmsh.model.occ.addLine(point_tags[i], point_tags[(i + 1) % 4]))

    loop = gmsh.model.occ.addCurveLoop(line_tags)
    surf = gmsh.model.occ.addPlaneSurface([loop])
    ext = gmsh.model.occ.extrude([(2, surf)], 0.0, 0.0, h)
    vol_tags = [tag for dim, tag in ext if dim == 3]

    return vol_tags[0], point_tags


def _collect_periodic_surfaces(final_vols: list[tuple[int, int]], a: float, tol: float = 2e-2) -> dict[str, list[int]]:
    side_groups = {"f1m": [], "f1p": [], "f2m": [], "f2p": []}

    boundaries = gmsh.model.getBoundary(final_vols, oriented=False, recursive=False)
    for dim, tag in boundaries:
        if dim != 2:
            continue

        cx, cy, _cz = gmsh.model.occ.getCenterOfMass(dim, tag)
        st = barycentric_coords(np.array([[cx, cy]], dtype=float), a=a)[0]

        if abs(st[0]) < tol:
            side_groups["f1m"].append(tag)
        if abs(st[0] - 1.0) < tol:
            side_groups["f1p"].append(tag)
        if abs(st[1]) < tol:
            side_groups["f2m"].append(tag)
        if abs(st[1] - 1.0) < tol:
            side_groups["f2p"].append(tag)

    return {k: sorted(set(v)) for k, v in side_groups.items()}


def generate_unit_cell_mesh(
    out_path: Path,
    a: float,
    h: float,
    w_intra: float,
    w_inter: float,
    ratio: float,
    mesh_size: float,
) -> tuple[Path, Path]:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size * 0.5)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

        gmsh.model.add(f"unit_cell_r{ratio:.2f}")

        a1, a2 = lattice_vectors(a)
        A = np.array([0.0, 0.0], dtype=float)
        B = 0.5 * a1
        C = 0.5 * a2

        w_inter_eff = w_intra * ratio if ratio > 0.0 else w_inter

        cell_vol = _add_cell_prism(a=a, h=h, lc=mesh_size)

        beam_specs: list[tuple[np.ndarray, np.ndarray, float]] = [
            (A, B, w_intra),
            (B, C, w_intra),
            (C, A, w_intra),
            (A, A + 0.5 * a1, w_inter_eff),
            (B, B + 0.5 * a1, w_inter_eff),
            (B, B + 0.5 * (a2 - a1), w_inter_eff),
            (C, C + 0.5 * a2, w_inter_eff),
        ]

        beam_vols: list[tuple[int, int]] = []
        junction_points: list[int] = []
        for start, end, width in beam_specs:
            vol, pts = _add_beam_volume(start=start, end=end, width=width, h=h, lc=mesh_size * 0.5)
            beam_vols.append((3, vol))
            junction_points.extend(pts)

        fused = beam_vols
        while len(fused) > 1:
            left = fused[0]
            right = fused[1]
            out, _ = gmsh.model.occ.fuse([left], [right], removeObject=True, removeTool=True)
            fused = [d for d in out if d[0] == 3] + fused[2:]

        clipped, _ = gmsh.model.occ.intersect(fused, [(3, cell_vol)], removeObject=True, removeTool=False)
        final_vols = [d for d in clipped if d[0] == 3]

        gmsh.model.occ.synchronize()

        if junction_points:
            gmsh.model.mesh.setSize([(0, p) for p in sorted(set(junction_points))], mesh_size * 0.5)

        gmsh.model.addPhysicalGroup(3, [tag for _, tag in final_vols], tag=101)
        gmsh.model.setPhysicalName(3, 101, "solid")

        side_groups = _collect_periodic_surfaces(final_vols=final_vols, a=a)
        tag_map = {"f1m": 201, "f1p": 202, "f2m": 203, "f2p": 204}
        for name, surfs in side_groups.items():
            if surfs:
                gmsh.model.addPhysicalGroup(2, surfs, tag=tag_map[name])
                gmsh.model.setPhysicalName(2, tag_map[name], name)

        gmsh.model.mesh.generate(3)

        out_path.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(out_path))

    finally:
        gmsh.finalize()

    msh = meshio.read(out_path)
    xdmf_path = out_path.with_suffix(".xdmf")
    meshio.write(xdmf_path, msh)
    return out_path, xdmf_path


def _parse_args() -> argparse.Namespace:
    defaults = LatticeParams()

    parser = argparse.ArgumentParser(description="Generate a breathing-kagome unit-cell 3D mesh")
    parser.add_argument("--a", type=float, default=defaults.a)
    parser.add_argument("--h", type=float, default=defaults.h)
    parser.add_argument("--w-intra", type=float, default=defaults.w_intra)
    parser.add_argument("--w-inter", type=float, default=defaults.w_inter)
    parser.add_argument("--ratio", type=float, default=defaults.w_inter / defaults.w_intra)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--mesh-size", type=float, default=0.0005)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    out = args.out if args.out.suffix == ".msh" else args.out.with_suffix(".msh")
    msh_path, xdmf_path = generate_unit_cell_mesh(
        out_path=out,
        a=args.a,
        h=args.h,
        w_intra=args.w_intra,
        w_inter=args.w_inter,
        ratio=args.ratio,
        mesh_size=args.mesh_size,
    )
    print(f"Wrote: {msh_path}")
    print(f"Wrote: {xdmf_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
