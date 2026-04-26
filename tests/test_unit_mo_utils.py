"""Unit tests for dft_loc.utils.mo_utils — gap coverage."""

from __future__ import annotations

import pytest

from dft_loc.utils.mo_utils import (
    build_atom_mo_coefficients,
    get_top_atom_coeff,
    parse_mo_coefficients,
)


# ---------------------------------------------------------------------------
# parse_mo_coefficients()
# ---------------------------------------------------------------------------

class TestParseMoCoefficients:
    def test_all_virtual_mos_returns_empty(self):
        lines = [
            "MO # 1 Energy: 0.1 Occ: 0.0",
            "0 C 2pz  0.5",
        ]
        result = parse_mo_coefficients(lines)
        assert result == {}

    def test_custom_occupied_threshold(self):
        lines = [
            "MO # 1 Energy: -0.5 Occ: 1.0",  # below default 1.5
            "0 N 2pz  0.4",
            "MO # 2 Energy: -0.3 Occ: 2.0",
            "0 N 2pz  0.6",
        ]
        # With threshold=0.5, MO 1 (occ=1.0) should be included
        result = parse_mo_coefficients(lines, occupied_threshold=0.5)
        assert 1 in result
        assert 2 in result

    def test_lines_before_first_header_are_ignored(self):
        lines = [
            "Some header text",
            "0 C 2pz  0.9",  # no MO header yet — should be skipped
            "MO # 1 Energy: -1.0 Occ: 2.0",
            "0 C 2pz  0.3",
        ]
        result = parse_mo_coefficients(lines)
        assert result == {1: {"0-C": 0.3}}

    def test_coefficients_accumulated_as_absolute_values(self):
        lines = [
            "MO # 1 Energy: -1.0 Occ: 2.0",
            "0 C 2pz  -0.30",
            "0 C 2px   0.20",  # same atom, different orbital → accumulated
        ]
        result = parse_mo_coefficients(lines)
        assert result[1]["0-C"] == pytest.approx(0.50)

    def test_short_lines_are_skipped(self):
        lines = [
            "MO # 1 Energy: -1.0 Occ: 2.0",
            "0 C",  # only 2 parts — too short
            "0 O 2pz  0.5",
        ]
        result = parse_mo_coefficients(lines)
        assert "0-C" not in result.get(1, {})
        assert result[1]["0-O"] == pytest.approx(0.5)

    def test_non_integer_first_part_skipped(self):
        lines = [
            "MO # 1 Energy: -1.0 Occ: 2.0",
            "X C 2pz  0.9",  # 'X' is not a digit → skip
            "1 N 2pz  0.4",
        ]
        result = parse_mo_coefficients(lines)
        assert "X-C" not in result.get(1, {})
        assert result[1]["1-N"] == pytest.approx(0.4)


# ---------------------------------------------------------------------------
# build_atom_mo_coefficients()
# ---------------------------------------------------------------------------

class TestBuildAtomMoCoefficients:
    def test_sorted_by_coefficient_descending(self):
        mo_dict = {1: {"0-C": 0.3, "1-O": 0.8}, 2: {"0-C": 0.9}}
        result = build_atom_mo_coefficients(mo_dict)
        # 0-C entries: (2, 0.9) > (1, 0.3)
        assert result["0-C"][0] == (2, 0.9)
        assert result["0-C"][1] == (1, 0.3)

    def test_max_occupied_mos_trims_oldest_mos(self):
        mo_dict = {1: {"0-C": 0.5}, 2: {"0-C": 0.4}, 3: {"0-C": 0.6}}
        result = build_atom_mo_coefficients(mo_dict, max_occupied_mos=2)
        # Should only use MOs 2 and 3 (the last 2)
        mo_indices = {mo for mo, _ in result["0-C"]}
        assert mo_indices == {2, 3}

    def test_max_occupied_mos_zero_uses_all_mos(self):
        # Python: -0 == 0, so mo_levels[-0:] == mo_levels[0:] (all elements).
        # max_occupied_mos=0 is a no-op in the current implementation.
        mo_dict = {1: {"0-C": 0.5}}
        result_zero = build_atom_mo_coefficients(mo_dict, max_occupied_mos=0)
        result_none = build_atom_mo_coefficients(mo_dict, max_occupied_mos=None)
        assert result_zero == result_none

    def test_multiple_atoms(self):
        mo_dict = {1: {"0-C": 0.5, "1-N": 0.3}, 2: {"1-N": 0.7}}
        result = build_atom_mo_coefficients(mo_dict)
        assert "0-C" in result
        assert "1-N" in result
        # 1-N: (2, 0.7) > (1, 0.3)
        assert result["1-N"][0] == (2, 0.7)


# ---------------------------------------------------------------------------
# get_top_atom_coeff()
# ---------------------------------------------------------------------------

class TestGetTopAtomCoeff:
    @pytest.fixture()
    def atom_coeff(self):
        return {
            "0-C":  [(1, 0.9), (2, 0.8), (3, 0.7), (4, 0.6), (5, 0.5)],
            "1-H":  [(1, 0.7), (2, 0.6)],
            "2-B":  [(1, 0.9), (2, 0.8), (3, 0.7), (4, 0.6)],
            "3-Al": [(1, 0.8), (2, 0.7), (3, 0.6), (4, 0.5)],
            "4-N":  [(1, 0.9), (2, 0.8), (3, 0.7), (4, 0.6), (5, 0.5)],
        }

    def test_hydrogen_gets_one_coefficient(self, atom_coeff):
        result = get_top_atom_coeff(atom_coeff)
        assert len(result["1-H"]) == 1

    def test_boron_gets_three_coefficients(self, atom_coeff):
        result = get_top_atom_coeff(atom_coeff)
        assert len(result["2-B"]) == 3

    def test_aluminum_gets_three_coefficients(self, atom_coeff):
        result = get_top_atom_coeff(atom_coeff)
        assert len(result["3-Al"]) == 3

    def test_default_non_special_gets_four(self, atom_coeff):
        result = get_top_atom_coeff(atom_coeff)
        assert len(result["4-N"]) == 4

    def test_central_atom_gets_num_coeff(self, atom_coeff):
        result = get_top_atom_coeff(atom_coeff, central_atom="0-C", num_coeff=2)
        assert len(result["0-C"]) == 2

    def test_central_atom_minimum_one(self, atom_coeff):
        # Even if num_coeff=0, central atom should get at least 1
        result = get_top_atom_coeff(atom_coeff, central_atom="0-C", num_coeff=0)
        assert len(result["0-C"]) >= 1

    def test_no_central_atom_default_four(self, atom_coeff):
        result = get_top_atom_coeff(atom_coeff)
        # Carbon (non-special, no central_atom set) → 4 coefficients
        assert len(result["0-C"]) == 4
