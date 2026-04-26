"""Integration tests for DFT-LOC — multi-module pipelines.

These tests exercise the handoffs between modules without relying on real
example files so they remain fast and self-contained.
"""

from __future__ import annotations

import pytest

from dft_loc.analysis.bond_analysis import (
    classify_bonds,
    infer_bond_pair_indices,
    valency_atom,
)
from dft_loc.utils.geometry_utils import bond_distances, build_connectivity
from dft_loc.utils.io_utils import AtomRecord, Molecule, atom_labels, open_xyz
from dft_loc.utils.mo_utils import (
    build_atom_mo_coefficients,
    get_top_atom_coeff,
    parse_mo_coefficients,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_mol(*atoms: tuple[str, float, float, float]) -> Molecule:
    records = tuple(AtomRecord(s, x, y, z) for s, x, y, z in atoms)
    return Molecule(atom_count=len(records), title="test", atoms=records)


WATER_LOG_LINES = [
    "MO # 1 Energy: -1.35 Occ: 2.0",
    "0 O 2s   0.80",
    "0 O 2pz  0.50",
    "1 H 1s   0.30",
    "MO # 2 Energy: -0.73 Occ: 2.0",
    "0 O 2pz  0.60",
    "1 H 1s   0.70",
    "2 H 1s   0.70",
    "MO # 3 Energy: -0.60 Occ: 2.0",
    "0 O 2pz  0.40",
    "2 H 1s   0.50",
    "MO # 4 Energy:  0.20 Occ: 0.0",   # virtual — should be excluded
    "0 O 2pz  9.99",
]


# ---------------------------------------------------------------------------
# Integration: geometry pipeline
# ---------------------------------------------------------------------------

class TestGeometryPipeline:
    """build_connectivity → bond_distances stay consistent."""

    def test_water_connectivity_then_distances(self):
        mol = make_mol(
            ("O", 0.0, 0.0, 0.000),
            ("H", 0.0, 0.0, 0.960),
            ("H", 0.930, 0.0, -0.240),
        )
        conn = build_connectivity(mol)
        dists = bond_distances(mol, conn)

        # Both O-H bonds should be present
        assert any("0-O" in k and "1-H" in k for k in dists)
        assert any("0-O" in k and "2-H" in k for k in dists)

        # Distance values should be positive and physically reasonable
        for d in dists.values():
            assert 0.5 < d < 2.0

    def test_isolated_noble_gas_absent_from_distances(self):
        mol = make_mol(
            ("O", 0.0, 0.0, 0.0),
            ("H", 0.0, 0.0, 0.96),
            ("He", 100.0, 0.0, 0.0),
        )
        conn = build_connectivity(mol)
        dists = bond_distances(mol, conn)
        assert all("2-He" not in k for k in dists)

    def test_atom_labels_match_connectivity_keys(self):
        mol = make_mol(
            ("C", 0.0, 0.0, 0.0),
            ("N", 1.2, 0.0, 0.0),
        )
        labels = set(atom_labels(mol))
        conn = build_connectivity(mol)
        assert set(conn.keys()) == labels


# ---------------------------------------------------------------------------
# Integration: MO pipeline
# ---------------------------------------------------------------------------

class TestMoPipeline:
    """parse_mo_coefficients → build_atom_mo_coefficients → get_top_atom_coeff."""

    def test_occupancy_filter_propagates(self):
        mo_dict = parse_mo_coefficients(WATER_LOG_LINES)
        # MO 4 is virtual — should not appear
        assert 4 not in mo_dict
        assert set(mo_dict.keys()) == {1, 2, 3}

    def test_atom_coeff_keys_are_consistent_with_parsed_labels(self):
        mo_dict = parse_mo_coefficients(WATER_LOG_LINES)
        atom_coeff = build_atom_mo_coefficients(mo_dict)
        for label in atom_coeff:
            idx, sym = label.split("-", 1)
            assert idx.isdigit()
            assert sym.isalpha()

    def test_top_coeff_hydrogen_trimmed_to_one(self):
        mo_dict = parse_mo_coefficients(WATER_LOG_LINES)
        atom_coeff = build_atom_mo_coefficients(mo_dict)
        top = get_top_atom_coeff(atom_coeff)
        for label, pairs in top.items():
            if label.endswith("-H"):
                assert len(pairs) == 1, f"{label} should have 1 coefficient"

    def test_max_occupied_mos_limits_search_space(self):
        mo_dict = parse_mo_coefficients(WATER_LOG_LINES)
        full = build_atom_mo_coefficients(mo_dict)
        limited = build_atom_mo_coefficients(mo_dict, max_occupied_mos=1)

        # Limited version should have fewer or equal MO entries per atom
        for label in limited:
            assert len(limited[label]) <= len(full[label])


# ---------------------------------------------------------------------------
# Integration: geometry + MO → bond analysis
# ---------------------------------------------------------------------------

class TestGeometryMoBondPipeline:
    """Full DFT-LOC pipeline without touching example files."""

    @pytest.fixture()
    def water_mol(self):
        return make_mol(
            ("O", 0.0, 0.0, 0.000),
            ("H", 0.0, 0.0, 0.960),
            ("H", 0.930, 0.0, -0.240),
        )

    @pytest.fixture()
    def water_bond_types(self, water_mol):
        mo_dict = parse_mo_coefficients(WATER_LOG_LINES)
        atom_coeff = build_atom_mo_coefficients(mo_dict)
        top_coeff = get_top_atom_coeff(atom_coeff)
        conn = build_connectivity(water_mol)
        bond_indices = infer_bond_pair_indices(top_coeff, conn)
        return classify_bonds(bond_indices)

    def test_water_has_oh_bonds(self, water_bond_types):
        # Both O-H bonds should appear
        bond_keys = set(water_bond_types.keys())
        assert any("0-O" in k for k in bond_keys)

    def test_water_bond_types_are_valid(self, water_bond_types):
        for btype in water_bond_types.values():
            assert btype in {"single", "double", "triple", "unknown"}

    def test_valency_oxygen_at_least_two(self, water_bond_types):
        v = valency_atom("0-O", water_bond_types)
        assert v >= 2  # oxygen forms 2 O-H bonds

    def test_valency_hydrogen_at_least_one(self, water_bond_types):
        for label in ("1-H", "2-H"):
            assert valency_atom(label, water_bond_types) >= 1


# ---------------------------------------------------------------------------
# Integration: CLI argument parsing
# ---------------------------------------------------------------------------

class TestCliIntegration:
    """Tests that focus on the argument-parsing layer of the CLI."""

    def test_analyze_requires_xyz(self):
        from dft_loc.cli import _build_parser
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["analyze", "--log", "x.log"])

    def test_analyze_requires_log(self):
        from dft_loc.cli import _build_parser
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["analyze", "--xyz", "x.xyz"])

    def test_no_subcommand_none_command(self):
        from dft_loc.cli import _build_parser
        parser = _build_parser()
        args = parser.parse_args([])
        assert args.command is None

    def test_max_occupied_mos_parsed_as_int(self):
        from dft_loc.cli import _build_parser
        parser = _build_parser()
        args = parser.parse_args(["analyze", "--xyz", "a.xyz", "--log", "b.log", "--max-occupied-mos", "5"])
        assert args.max_occupied_mos == 5

    def test_analyze_with_tmp_files(self, tmp_path):
        """End-to-end: call main() with in-process tmp files."""
        from dft_loc.cli import main

        xyz = tmp_path / "mol.xyz"
        log = tmp_path / "mol.log"
        xyz.write_text(
            "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n",
            encoding="utf-8",
        )
        log.write_text(
            "\n".join(WATER_LOG_LINES),
            encoding="utf-8",
        )
        # Should not raise
        main(["analyze", "--xyz", str(xyz), "--log", str(log)])

    def test_analyze_with_max_occupied_mos_flag(self, tmp_path, capsys):
        from dft_loc.cli import main

        xyz = tmp_path / "mol.xyz"
        log = tmp_path / "mol.log"
        xyz.write_text(
            "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n",
            encoding="utf-8",
        )
        log.write_text("\n".join(WATER_LOG_LINES), encoding="utf-8")

        main(["analyze", "--xyz", str(xyz), "--log", str(log), "--max-occupied-mos", "1"])
        captured = capsys.readouterr()
        assert "Bond summary" in captured.out
