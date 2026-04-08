from __future__ import annotations

from pathlib import Path

import pytest

from dft_loc.analysis.bond_analysis import bond_order1, classify_bonds, infer_bond_pair_indices, valency_atom
from dft_loc.utils.geometry_utils import angle, bond_distances, build_connectivity, distance
from dft_loc.utils.io_utils import AtomRecord, Molecule, atom_labels, open_xyz
from dft_loc.utils.mo_utils import build_atom_mo_coefficients, get_top_atom_coeff, parse_mo_coefficients


def make_molecule(*atoms: tuple[str, float, float, float], title: str = "test") -> Molecule:
    records = tuple(AtomRecord(symbol, x, y, z) for symbol, x, y, z in atoms)
    return Molecule(atom_count=len(records), title=title, atoms=records)


def test_open_xyz_skips_charge_multiplicity_line(tmp_path: Path):
    xyz_path = tmp_path / "water.xyz"
    xyz_path.write_text(
        "\n".join(
            [
                "3",
                "water",
                "0 1",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 0.9",
                "H 0.0 0.8 0.0",
            ]
        ),
        encoding="utf-8",
    )

    mol = open_xyz(xyz_path)

    assert mol.atom_count == 3
    assert atom_labels(mol) == ("0-O", "1-H", "2-H")


def test_open_xyz_rejects_wrong_atom_count(tmp_path: Path):
    xyz_path = tmp_path / "bad.xyz"
    xyz_path.write_text(
        "\n".join(
            [
                "2",
                "broken",
                "H 0.0 0.0 0.0",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="expected 2 atom lines, found 1"):
        open_xyz(xyz_path)


def test_geometry_helpers_build_connectivity_and_distances():
    mol = make_molecule(
        ("O", 0.0, 0.0, 0.0),
        ("H", 0.0, 0.0, 0.96),
        ("H", 0.93, 0.0, -0.24),
        ("He", 5.0, 5.0, 5.0),
    )

    connectivity = build_connectivity(mol)
    distances = bond_distances(mol, connectivity)

    assert connectivity["0-O"] == ["1-H", "2-H"]
    assert connectivity["3-He"] == []
    assert pytest.approx(distance(mol, 0, 1), rel=1e-6) == 0.96
    assert pytest.approx(angle(mol, 1, 0, 2), rel=1e-6) == 104.4702941000659
    assert "0-O:1-H" in distances
    assert "0-O:2-H" in distances
    assert all("3-He" not in pair for pair in distances)


def test_angle_raises_for_zero_length_vector():
    mol = make_molecule(
        ("H", 0.0, 0.0, 0.0),
        ("H", 0.0, 0.0, 0.0),
        ("H", 1.0, 0.0, 0.0),
    )

    with pytest.raises(ValueError, match="zero-length vector"):
        angle(mol, 0, 1, 2)


def test_parse_mo_coefficients_sums_absolute_values_and_filters_virtual_mos():
    lines = [
        "MO # 1 Energy: -0.8 Occ: 2.0",
        "0 C 2pz    -0.30",
        "0 C 2px     0.20",
        "1 O 2pz     0.40",
        "MO # 2 Energy: -0.1 Occ: 0.0",
        "0 C 2pz     9.99",
    ]

    mo_dict = parse_mo_coefficients(lines)

    assert mo_dict == {1: {"0-C": 0.5, "1-O": 0.4}}


def test_build_atom_mo_coefficients_and_top_coefficients():
    mo_dict = {
        1: {"0-C": 0.1, "1-H": 0.7},
        2: {"0-C": 0.9, "1-H": 0.2, "2-B": 0.8},
        3: {"0-C": 0.4, "2-B": 0.6},
        4: {"0-C": 0.3},
        5: {"0-C": 0.2},
    }

    atom_coeff = build_atom_mo_coefficients(mo_dict, max_occupied_mos=4)
    top_coeff = get_top_atom_coeff(atom_coeff, central_atom="0-C", num_coeff=2)

    assert atom_coeff["0-C"] == [(2, 0.9), (3, 0.4), (4, 0.3), (5, 0.2)]
    assert top_coeff["0-C"] == [(2, 0.9), (3, 0.4)]
    assert top_coeff["1-H"] == [(2, 0.2)]
    assert top_coeff["2-B"] == [(2, 0.8), (3, 0.6)]


def test_bond_analysis_helpers_cover_multiple_orders():
    top_atom_coeff = {
        "0-C": [(1, 0.9), (2, 0.8), (3, 0.7), (4, 0.6)],
        "1-O": [(2, 0.9), (3, 0.8)],
        "2-N": [(1, 0.9), (2, 0.7), (3, 0.6)],
        "3-H": [(4, 0.5)],
    }
    connectivity = {
        "0-C": ["1-O", "2-N", "3-H"],
        "1-O": ["0-C"],
        "2-N": ["0-C"],
        "3-H": ["0-C"],
    }

    indices = infer_bond_pair_indices(top_atom_coeff, connectivity)
    bond_types = classify_bonds(indices)

    assert indices["0-C:1-O"] == [2, 3]
    assert indices["0-C:2-N"] == [1, 2, 3]
    assert indices["0-C:3-H"] == [4]
    assert bond_types == {
        "0-C:1-O": "double",
        "0-C:2-N": "triple",
        "0-C:3-H": "single",
    }
    assert valency_atom("0-C", bond_types) == 6
    assert bond_order1("0-C", bond_types) == (1, 1, 1)
