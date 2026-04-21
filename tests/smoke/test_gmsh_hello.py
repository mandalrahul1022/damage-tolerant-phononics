"""
Smoke test: verify the Gmsh Python API can build and mesh a 2D square.

Implements: BLUEPRINT.md Section 5.1 (Phase 2 Step 5.1, Step A5 of the
session plan). Confirms the Gmsh binding is functional before we use it
in Phase 3 to build the parametric 3D kagome unit-cell mesh.

Verifies:
    Creates a unit square via the OpenCASCADE kernel, meshes it with
    characteristic length 0.1, and counts the resulting 2D elements
    (triangles). Passes if the element count is strictly positive.
    Typical count at cl=0.1 is ~200 triangles.

How to run:
    conda activate kagome-hoti
    python tests/smoke/test_gmsh_hello.py
or:
    pytest tests/smoke/test_gmsh_hello.py -v

Expected output on PASS:
    "Gmsh 2D triangle count: NNN"
    "PASS". Exit code 0.

Pass/fail criteria:
    PASS: gmsh initializes, meshes the square, reports >0 triangles, exit 0.
    FAIL: any exception, zero elements, or nonzero exit.
"""

from __future__ import annotations

import sys

import gmsh


def _count_square_triangles() -> int:
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("unit_square")
        gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, 1.0, 1.0)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1)
        gmsh.model.mesh.generate(2)

        element_types, element_tags, _ = gmsh.model.mesh.getElements(dim=2)
        return sum(len(tags) for tags in element_tags)
    finally:
        gmsh.finalize()


def test_square_has_triangles() -> None:
    count = _count_square_triangles()
    assert count > 0, f"Gmsh produced {count} 2D elements; expected > 0"


if __name__ == "__main__":
    count = _count_square_triangles()
    print(f"Gmsh 2D triangle count: {count}")
    if count > 0:
        print("PASS")
        sys.exit(0)
    print("FAIL")
    sys.exit(1)
