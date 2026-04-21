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


def _build_bond_specs(a: float, w_intra: float, w_inter: float) -> list[tuple[np.ndarray, np.ndarray, float]]:
    """Build a deduplicated set of breathing-kagome bonds around the unit cell.

    We include a 3x3 neighborhood of lattice translations and clip against the
    unit-cell prism. This guarantees periodic boundary cuts on opposite sides.
    """
    a1, a2 = lattice_vectors(a)

    def _sites(R: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        A = R
        B = R + 0.5 * a1
        C = R + 0.5 * a2
        return A, B, C

    def _bond_key(p: np.ndarray, q: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
        p_key = tuple(np.round(p, 10))
        q_key = tuple(np.round(q, 10))
        return tuple(sorted((p_key, q_key)))

    specs: list[tuple[np.ndarray, np.ndarray, float]] = []
    seen: set[tuple[tuple[float, float], tuple[float, float]]] = set()

    for i in (-1, 0, 1):
        for j in (0, 1):
            R = i * a1 + j * a2
            A, B, C = _sites(R)

            bonds = [
                (A, B, w_intra),
                (B, C, w_intra),
                (C, A, w_intra),
                (A, B + a1, w_inter),
                (B, C + (a2 - a1), w_inter),
                (A, C + a2, w_inter),
            ]

            for p, q, w in bonds:
                k = _bond_key(p, q)
                if k in seen:
                    continue
                seen.add(k)
                specs.append((np.asarray(p, dtype=float), np.asarray(q, dtype=float), float(w)))

    return specs


def _collect_periodic_surfaces(final_vols: list[tuple[int, int]], a: float, h: float, tol: float = 2e-2) -> dict[str, list[int]]:
    side_groups: dict[str, list[int]] = {"f1m": [], "f1p": [], "f2m": [], "f2p": []}
    plane_tol = 1e-8

    boundaries = gmsh.model.getBoundary(final_vols, oriented=False, recursive=False)
    for dim, tag in boundaries:
        if dim != 2:
            continue

        _xmin, _ymin, zmin, _xmax, _ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        if (zmax - zmin) < 0.9 * h:
            continue

        # Use all surface vertices: periodic faces are clipping faces on unit-cell
        # boundary planes s=0/1 and t=0/1 in lattice coordinates.
        pts = []
        sub = gmsh.model.getBoundary([(2, tag)], oriented=False, recursive=True)
        for dsub, tsub in sub:
            if dsub != 0:
                continue
            x, y, _z = gmsh.model.getValue(0, tsub, [])
            pts.append([x, y])
        if not pts:
            continue
        st = barycentric_coords(np.asarray(pts, dtype=float), a=a)
        s = st[:, 0]
        t = st[:, 1]
        dists = {
            "f1m": float(np.max(np.abs(s))),
            "f1p": float(np.max(np.abs(s - 1.0))),
            "f2m": float(np.max(np.abs(t))),
            "f2p": float(np.max(np.abs(t - 1.0))),
        }
        ordered = sorted(dists.items(), key=lambda kv: kv[1])
        (best_side, best_d), (_second_side, second_d) = ordered[0], ordered[1]

        if best_d > plane_tol:
            continue

        if second_d < 10.0 * best_d:
            cx, cy, _cz = gmsh.model.occ.getCenterOfMass(dim, tag)
            raise RuntimeError(
                f"Corner ambiguity for face tag={tag} at centroid=({cx:+.6e},{cy:+.6e}): "
                f"distances f1m={dists['f1m']:.3e} f1p={dists['f1p']:.3e} "
                f"f2m={dists['f2m']:.3e} f2p={dists['f2p']:.3e}"
            )

        side_groups[best_side].append(tag)

    _ = tol  # retained for signature stability; classification is tolerance-free
    side_groups = {k: sorted(set(v)) for k, v in side_groups.items()}
    for side in ("f1m", "f1p", "f2m", "f2p"):
        if len(side_groups[side]) != 1:
            raise RuntimeError(
                f"Periodic side {side} expected 1 face, found {len(side_groups[side])}: {side_groups[side]}"
            )
    return side_groups


def _apply_periodic(master_tags: list[int], slave_tags: list[int], translation_xy: np.ndarray, side_name: str) -> None:
    """Pair each slave face to its master by centroid-after-translation, then call setPeriodic per pair."""
    if len(master_tags) != len(slave_tags):
        raise RuntimeError(
            f"Periodic side count mismatch for {side_name}: "
            f"|masters|={len(master_tags)} (tags={master_tags}), "
            f"|slaves|={len(slave_tags)} (tags={slave_tags})"
        )

    tx, ty = float(translation_xy[0]), float(translation_xy[1])
    affine = [
        1.0, 0.0, 0.0, tx,
        0.0, 1.0, 0.0, ty,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ]

    master_centroids = {
        t: np.array(gmsh.model.occ.getCenterOfMass(2, t)[:2], dtype=float)
        for t in master_tags
    }

    matched: set[int] = set()
    for slave_tag in slave_tags:
        sc = np.array(gmsh.model.occ.getCenterOfMass(2, slave_tag)[:2], dtype=float)
        expected_master = sc - np.array([tx, ty], dtype=float)

        best_tag = None
        best_d = float("inf")
        for mt, mc in master_centroids.items():
            if mt in matched:
                continue
            d = float(np.linalg.norm(mc - expected_master))
            if d < best_d:
                best_d = d
                best_tag = mt

        if best_tag is None or best_d > 1e-5:
            raise RuntimeError(
                f"setPeriodic pairing failed for {side_name}: slave surface {slave_tag} "
                f"at centroid {sc.tolist()} has no master within 1e-5 m "
                f"(best master {best_tag}, best dist {best_d:.2e}, translation {translation_xy.tolist()})"
            )
        matched.add(best_tag)
        gmsh.model.mesh.setPeriodic(2, [slave_tag], [best_tag], affine)


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

        beam_specs = _build_bond_specs(a=a, w_intra=w_intra, w_inter=w_inter_eff)

        beam_vols: list[tuple[int, int]] = []
        junction_points: list[int] = []
        for start, end, width in beam_specs:
            vol, pts = _add_beam_volume(start=start, end=end, width=width, h=h, lc=mesh_size * 0.5)
            beam_vols.append((3, vol))
            junction_points.extend(pts)

        clipped, _ = gmsh.model.occ.intersect(beam_vols, [(3, cell_vol)], removeObject=True, removeTool=False)
        clipped_vols = [d for d in clipped if d[0] == 3]
        if not clipped_vols:
            raise RuntimeError("Boolean intersection produced no clipped beam volumes")
        if len(clipped_vols) == 1:
            final_vols = clipped_vols
        else:
            fused, _ = gmsh.model.occ.fuse([clipped_vols[0]], clipped_vols[1:], removeObject=True, removeTool=True)
            final_vols = [d for d in fused if d[0] == 3]

        gmsh.model.occ.synchronize()

        if junction_points:
            gmsh.model.mesh.setSize([(0, p) for p in sorted(set(junction_points))], mesh_size * 0.5)

        gmsh.model.addPhysicalGroup(3, [tag for _, tag in final_vols], tag=101)
        gmsh.model.setPhysicalName(3, 101, "solid")

        side_groups = _collect_periodic_surfaces(final_vols=final_vols, a=a, h=h)
        tag_map = {"f1m": 201, "f1p": 202, "f2m": 203, "f2p": 204}
        for name, surfs in side_groups.items():
            if surfs:
                gmsh.model.addPhysicalGroup(2, surfs, tag=tag_map[name])
                gmsh.model.setPhysicalName(2, tag_map[name], name)

        a1_vec, a2_vec = lattice_vectors(a)
        _apply_periodic(master_tags=side_groups["f1m"], slave_tags=side_groups["f1p"], translation_xy=a1_vec, side_name="f1m->f1p (+a1)")
        _apply_periodic(master_tags=side_groups["f2m"], slave_tags=side_groups["f2p"], translation_xy=a2_vec, side_name="f2m->f2p (+a2)")

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
    parser.add_argument("--mesh-size", type=float, default=0.00052)
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
