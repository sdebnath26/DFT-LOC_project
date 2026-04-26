"""Unit tests for dft_loc.utils.geometry_utils — gap coverage."""

from __future__ import annotations

import math

import pytest

from dft_loc.utils.geometry_utils import (
    COVALENT_RADII,
    angle,
    bond_distances,
    build_connectivity,
    distance,
)
from dft_loc.utils.io_utils import AtomRecord, Molecule


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_mol(*atoms: tuple[str, float, float, float]) -> Molecule:
    records = tuple(AtomRecord(s, x, y, z) for s, x, y, z in atoms)
    return Molecule(atom_count=len(records), title="test", atoms=records)


# ---------------------------------------------------------------------------
# distance()
# ---------------------------------------------------------------------------

class TestDistance:
    def test_same_atom_is_zero(self):
        mol = make_mol(("C", 1.0, 2.0, 3.0))
        assert distance(mol, 0, 0) == 0.0

    def test_unit_x_axis(self):
        mol = make_mol(("H", 0.0, 0.0, 0.0), ("H", 1.0, 0.0, 0.0))
        assert distance(mol, 0, 1) == pytest.approx(1.0)

    def test_3d_pythagoras(self):
        mol = make_mol(("C", 0.0, 0.0, 0.0), ("O", 1.0, 1.0, 1.0))
        assert distance(mol, 0, 1) == pytest.approx(math.sqrt(3), rel=1e-9)

    def test_symmetric(self):
        mol = make_mol(("N", 0.3, -0.7, 1.2), ("S", -1.1, 2.3, 0.0))
        assert distance(mol, 0, 1) == pytest.approx(distance(mol, 1, 0))


# ---------------------------------------------------------------------------
# angle()
# ---------------------------------------------------------------------------

class TestAngle:
    def test_right_angle(self):
        # H-O-H at 90°
        mol = make_mol(
            ("H", 1.0, 0.0, 0.0),  # 0 – endpoint a
            ("O", 0.0, 0.0, 0.0),  # 1 – central b
            ("H", 0.0, 1.0, 0.0),  # 2 – endpoint c
        )
        assert angle(mol, 0, 1, 2) == pytest.approx(90.0, rel=1e-9)

    def test_linear_molecule_180_degrees(self):
        # H-C-H collinear along x
        mol = make_mol(
            ("H", -1.0, 0.0, 0.0),
            ("C",  0.0, 0.0, 0.0),
            ("H",  1.0, 0.0, 0.0),
        )
        assert angle(mol, 0, 1, 2) == pytest.approx(180.0, rel=1e-9)

    def test_zero_length_vector_raises(self):
        mol = make_mol(
            ("H", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 0.0),  # coincident with atom 0
            ("H", 1.0, 0.0, 0.0),
        )
        with pytest.raises(ValueError, match="zero-length vector"):
            angle(mol, 0, 1, 2)

    def test_known_water_angle(self):
        # Same geometry used in test_unit_helpers.py
        mol = make_mol(
            ("O", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 0.96),
            ("H", 0.93, 0.0, -0.24),
        )
        assert angle(mol, 1, 0, 2) == pytest.approx(104.47, rel=1e-3)


# ---------------------------------------------------------------------------
# build_connectivity()
# ---------------------------------------------------------------------------

class TestBuildConnectivity:
    def test_isolated_atom_has_empty_neighbors(self):
        mol = make_mol(
            ("He", 0.0, 0.0, 0.0),
            ("Ne", 100.0, 0.0, 0.0),
        )
        conn = build_connectivity(mol)
        assert conn["0-He"] == []
        assert conn["1-Ne"] == []

    def test_all_atoms_present_even_if_isolated(self):
        mol = make_mol(("Ar", 0.0, 0.0, 0.0))
        conn = build_connectivity(mol)
        assert "0-Ar" in conn

    def test_custom_tolerance_tighter_breaks_bond(self):
        # O-H at 0.96 Å; with tight tolerance, covalent sum = 0.66+0.31 = 0.97
        # Default tol=0.40 → threshold 1.37 → bonded.
        # With tol=-0.05 → threshold 0.92 → not bonded.
        mol = make_mol(
            ("O", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 0.96),
        )
        tight = build_connectivity(mol, tolerance=-0.05)
        default = build_connectivity(mol)
        assert default["0-O"] == ["1-H"]
        assert tight["0-O"] == []

    def test_unknown_element_uses_fallback_radius(self):
        # Xe is not in COVALENT_RADII; should use 0.75 fallback
        assert "Xe" not in COVALENT_RADII
        mol = make_mol(
            ("Xe", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 1.0),
        )
        # 0.75 + 0.31 + 0.40 = 1.46 → 1.0 Å is within threshold → bonded
        conn = build_connectivity(mol)
        assert "1-H" in conn["0-Xe"]

    def test_bidirectional(self):
        mol = make_mol(
            ("C", 0.0, 0.0, 0.0),
            ("H", 1.0, 0.0, 0.0),
        )
        conn = build_connectivity(mol)
        assert "1-H" in conn["0-C"]
        assert "0-C" in conn["1-H"]


# ---------------------------------------------------------------------------
# bond_distances()
# ---------------------------------------------------------------------------

class TestBondDistances:
    def test_empty_connectivity_returns_empty(self):
        mol = make_mol(("C", 0.0, 0.0, 0.0))
        conn = {"0-C": []}
        assert bond_distances(mol, conn) == {}

    def test_each_pair_appears_once(self):
        mol = make_mol(
            ("O", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 0.96),
            ("H", 0.93, 0.0, -0.24),
        )
        conn = build_connectivity(mol)
        distances = bond_distances(mol, conn)
        # No duplicate keys
        assert len(distances) == len(set(distances))

    def test_distance_value_is_correct(self):
        mol = make_mol(
            ("C", 0.0, 0.0, 0.0),
            ("O", 0.0, 0.0, 1.2),
        )
        conn = {"0-C": ["1-O"], "1-O": ["0-C"]}
        dists = bond_distances(mol, conn)
        assert len(dists) == 1
        assert list(dists.values())[0] == pytest.approx(1.2, rel=1e-9)
