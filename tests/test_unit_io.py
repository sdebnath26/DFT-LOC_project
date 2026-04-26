"""Unit tests for dft_loc.utils.io_utils — gap coverage."""

from __future__ import annotations

import pytest

from dft_loc.utils.io_utils import (
    AtomRecord,
    Molecule,
    atom_labels,
    element_symbols,
    open_xyz,
)


# ---------------------------------------------------------------------------
# open_xyz() edge cases
# ---------------------------------------------------------------------------

class TestOpenXyz:
    def test_file_too_short_raises(self, tmp_path):
        p = tmp_path / "short.xyz"
        p.write_text("1\n", encoding="utf-8")
        with pytest.raises(ValueError, match="too short"):
            open_xyz(p)

    def test_malformed_atom_row_raises(self, tmp_path):
        p = tmp_path / "bad_atom.xyz"
        p.write_text("1\ntitle\nH 0.0 0.0\n", encoding="utf-8")  # only 3 parts
        with pytest.raises(ValueError):
            open_xyz(p)

    def test_non_integer_atom_count_raises(self, tmp_path):
        p = tmp_path / "bad_count.xyz"
        p.write_text("two\ntitle\nH 0.0 0.0 0.0\n", encoding="utf-8")
        with pytest.raises((ValueError, TypeError)):
            open_xyz(p)

    def test_charge_multiplicity_line_stripped(self, tmp_path):
        p = tmp_path / "cm.xyz"
        p.write_text(
            "2\nwater\n0 1\nO 0.0 0.0 0.0\nH 0.0 0.0 0.9\n",
            encoding="utf-8",
        )
        mol = open_xyz(p)
        assert mol.atom_count == 2
        assert mol.atoms[0].symbol == "O"
        assert mol.atoms[1].symbol == "H"

    def test_correct_coordinates_parsed(self, tmp_path):
        p = tmp_path / "co.xyz"
        p.write_text("2\nCO\nC 0.0 0.0 0.0\nO 0.0 0.0 1.128\n", encoding="utf-8")
        mol = open_xyz(p)
        assert mol.atoms[0].symbol == "C"
        assert mol.atoms[1].z == pytest.approx(1.128)

    def test_title_preserved(self, tmp_path):
        p = tmp_path / "t.xyz"
        p.write_text("1\nmy molecule\nH 0 0 0\n", encoding="utf-8")
        mol = open_xyz(p)
        assert mol.title == "my molecule"

    def test_frozen_dataclass_immutable(self, tmp_path):
        p = tmp_path / "h.xyz"
        p.write_text("1\ntitle\nH 0.0 0.0 0.0\n", encoding="utf-8")
        mol = open_xyz(p)
        with pytest.raises((TypeError, AttributeError)):
            mol.atom_count = 99  # type: ignore[misc]


# ---------------------------------------------------------------------------
# atom_labels() / element_symbols()
# ---------------------------------------------------------------------------

class TestLabelsAndSymbols:
    @pytest.fixture()
    def simple_mol(self):
        return Molecule(
            atom_count=3,
            title="test",
            atoms=(
                AtomRecord("C", 0.0, 0.0, 0.0),
                AtomRecord("Cl", 1.0, 0.0, 0.0),
                AtomRecord("H", 2.0, 0.0, 0.0),
            ),
        )

    def test_atom_labels_format(self, simple_mol):
        labels = atom_labels(simple_mol)
        assert labels == ("0-C", "1-Cl", "2-H")

    def test_element_symbols(self, simple_mol):
        syms = element_symbols(simple_mol)
        assert syms == ("C", "Cl", "H")

    def test_single_atom(self):
        mol = Molecule(
            atom_count=1,
            title="t",
            atoms=(AtomRecord("N", 0.0, 0.0, 0.0),),
        )
        assert atom_labels(mol) == ("0-N",)
        assert element_symbols(mol) == ("N",)

    def test_empty_molecule(self):
        mol = Molecule(atom_count=0, title="empty", atoms=())
        assert atom_labels(mol) == ()
        assert element_symbols(mol) == ()
