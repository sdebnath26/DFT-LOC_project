from __future__ import annotations

from collections import Counter


def _common_mos(a: list[tuple[int, float]], b: list[tuple[int, float]]) -> list[int]:
    ma = {mo for mo, _ in a}
    mb = {mo for mo, _ in b}
    return sorted(ma.intersection(mb))


def infer_bond_pair_indices(
    top_atom_coeff: dict[str, list[tuple[int, float]]],
    connectivity: dict[str, list[str]],
) -> dict[str, list[int]]:
    """For each connected atom pair, compute intersecting MO indices."""
    out: dict[str, list[int]] = {}
    seen: set[tuple[str, str]] = set()

    for a, nbrs in connectivity.items():
        for b in nbrs:
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            out[f"{a}:{b}"] = _common_mos(top_atom_coeff.get(a, []), top_atom_coeff.get(b, []))

    return out


def classify_bonds(bond_pair_indices: dict[str, list[int]]) -> dict[str, str]:
    """Classify by count of shared MOs as single/double/triple/unknown."""
    out: dict[str, str] = {}
    for pair, indices in bond_pair_indices.items():
        n = len(indices)
        if n <= 1:
            out[pair] = "single"
        elif n == 2:
            out[pair] = "double"
        elif n == 3:
            out[pair] = "triple"
        else:
            out[pair] = "unknown"
    return out


def valency_atom(atom_label: str, bond_types: dict[str, str]) -> int:
    """Compute simple valency-like count from bond type weights."""
    order = {"single": 1, "double": 2, "triple": 3, "unknown": 1}
    total = 0
    for pair, typ in bond_types.items():
        a, b = pair.split(":")
        if atom_label == a or atom_label == b:
            total += order[typ]
    return total


def bond_order1(atom_label: str, bond_types: dict[str, str]) -> tuple[int, int, int]:
    """Legacy helper name retained: returns (single, double, triple) counts for an atom."""
    c = Counter()
    for pair, typ in bond_types.items():
        a, b = pair.split(":")
        if atom_label == a or atom_label == b:
            c[typ] += 1
    return c["single"], c["double"], c["triple"]
