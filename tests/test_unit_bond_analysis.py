"""Unit tests for dft_loc.analysis.bond_analysis — gap coverage."""

from __future__ import annotations

import pytest

from dft_loc.analysis.bond_analysis import (
    bond_order1,
    classify_bonds,
    infer_bond_pair_indices,
    valency_atom,
)


# ---------------------------------------------------------------------------
# classify_bonds()
# ---------------------------------------------------------------------------

class TestClassifyBonds:
    def test_zero_shared_mos_is_single(self):
        assert classify_bonds({"0-C:1-H": []}) == {"0-C:1-H": "single"}

    def test_one_shared_mo_is_single(self):
        assert classify_bonds({"0-C:1-H": [1]}) == {"0-C:1-H": "single"}

    def test_two_shared_mos_is_double(self):
        assert classify_bonds({"0-C:1-O": [1, 2]}) == {"0-C:1-O": "double"}

    def test_three_shared_mos_is_triple(self):
        assert classify_bonds({"0-C:1-N": [1, 2, 3]}) == {"0-C:1-N": "triple"}

    def test_four_or_more_is_unknown(self):
        result = classify_bonds({"0-C:1-X": [1, 2, 3, 4]})
        assert result == {"0-C:1-X": "unknown"}

    def test_empty_input(self):
        assert classify_bonds({}) == {}

    def test_mixed_bond_types(self):
        result = classify_bonds({
            "0-C:1-H": [5],
            "0-C:2-O": [1, 2],
            "0-C:3-N": [1, 2, 3],
        })
        assert result["0-C:1-H"] == "single"
        assert result["0-C:2-O"] == "double"
        assert result["0-C:3-N"] == "triple"


# ---------------------------------------------------------------------------
# valency_atom()
# ---------------------------------------------------------------------------

class TestValencyAtom:
    def test_no_bonds_returns_zero(self):
        bond_types = {"1-O:2-H": "single"}
        assert valency_atom("0-C", bond_types) == 0

    def test_single_bond(self):
        bond_types = {"0-C:1-H": "single"}
        assert valency_atom("0-C", bond_types) == 1

    def test_double_bond_counts_two(self):
        bond_types = {"0-C:1-O": "double"}
        assert valency_atom("0-C", bond_types) == 2

    def test_triple_bond_counts_three(self):
        bond_types = {"0-C:1-N": "triple"}
        assert valency_atom("0-C", bond_types) == 3

    def test_unknown_bond_counts_one(self):
        bond_types = {"0-C:1-X": "unknown"}
        assert valency_atom("0-C", bond_types) == 1

    def test_multiple_bonds_summed(self):
        bond_types = {
            "0-C:1-O": "double",   # 2
            "0-C:2-N": "triple",   # 3
            "0-C:3-H": "single",   # 1
        }
        assert valency_atom("0-C", bond_types) == 6

    def test_atom_as_second_in_pair(self):
        bond_types = {"0-C:1-O": "double"}
        assert valency_atom("1-O", bond_types) == 2


# ---------------------------------------------------------------------------
# bond_order1()
# ---------------------------------------------------------------------------

class TestBondOrder1:
    def test_no_bonds_returns_zeros(self):
        assert bond_order1("0-C", {}) == (0, 0, 0)

    def test_mixed_bonds_counted_correctly(self):
        bond_types = {
            "0-C:1-H": "single",
            "0-C:2-H": "single",
            "0-C:3-O": "double",
            "0-C:4-N": "triple",
        }
        singles, doubles, triples = bond_order1("0-C", bond_types)
        assert singles == 2
        assert doubles == 1
        assert triples == 1

    def test_atom_as_second_in_pair(self):
        bond_types = {"0-C:1-H": "single"}
        singles, doubles, triples = bond_order1("1-H", bond_types)
        assert singles == 1
        assert doubles == 0
        assert triples == 0


# ---------------------------------------------------------------------------
# infer_bond_pair_indices()
# ---------------------------------------------------------------------------

class TestInferBondPairIndices:
    def test_empty_connectivity_returns_empty(self):
        result = infer_bond_pair_indices({}, {})
        assert result == {}

    def test_isolated_atoms_produce_no_pairs(self):
        top = {"0-C": [(1, 0.9)], "1-O": [(2, 0.8)]}
        conn = {"0-C": [], "1-O": []}
        result = infer_bond_pair_indices(top, conn)
        assert result == {}

    def test_missing_atom_in_top_coeff_yields_empty_intersection(self):
        # Atom in connectivity but not in top_atom_coeff
        top = {"0-C": [(1, 0.9), (2, 0.8)]}
        conn = {"0-C": ["1-H"], "1-H": ["0-C"]}
        result = infer_bond_pair_indices(top, conn)
        assert result["0-C:1-H"] == []

    def test_each_pair_appears_once(self):
        top = {"0-C": [(1, 0.9)], "1-O": [(1, 0.8)]}
        conn = {"0-C": ["1-O"], "1-O": ["0-C"]}
        result = infer_bond_pair_indices(top, conn)
        assert len(result) == 1

    def test_shared_mo_detected(self):
        top = {
            "0-C": [(1, 0.9), (2, 0.8)],
            "1-O": [(2, 0.7), (3, 0.6)],
        }
        conn = {"0-C": ["1-O"], "1-O": ["0-C"]}
        result = infer_bond_pair_indices(top, conn)
        assert result["0-C:1-O"] == [2]  # only MO 2 is shared
