from __future__ import annotations

from analysis.loc_corrections import (
    ATOM_PARAMS,
    BOND_PARAMS,
    ENV_PARAMS,
    LocAtom,
    LocBond,
    LocCorrectionBreakdown,
)


def atomic_loc_correction(atom: LocAtom) -> float:
    e = atom.element
    h = atom.hybridization

    if atom.octet_expanded and e in {"Cl", "P", "S"}:
        return ATOM_PARAMS["OCT_EXP"]

    if e == "Be" and h == "sp":
        return ATOM_PARAMS["Be_sp"]

    if e in {"N", "P"}:
        if h == "sp":
            return ATOM_PARAMS[f"{e}_sp"]
        if h == "sp2":
            return ATOM_PARAMS[f"{e}_sp2"]
        if h == "sp3":
            return ATOM_PARAMS[f"{e}_sp3"]

    if e == "O":
        if h == "sp2":
            return ATOM_PARAMS["O_sp2"]
        if h == "sp3":
            return ATOM_PARAMS["O_sp3"]

    return 0.0


def bond_loc_correction(bond: LocBond) -> float:
    if bond.charge_transfer:
        return BOND_PARAMS["CT"]
    if bond.polarized_class is not None:
        return BOND_PARAMS[bond.polarized_class]
    if bond.multiple_class is not None:
        return BOND_PARAMS[bond.multiple_class]
    if bond.length_class is not None:
        return BOND_PARAMS[bond.length_class]
    return 0.0


def environment_loc_correction(
    bond_key: str,
    bond: LocBond,
    bond_map: dict[str, LocBond],
) -> float:
    # LOC ESBC applies to single bonds between non-H/F atoms,
    # excluding special cases like small rings.
    ae = bond.a.split("-")[1]
    be = bond.b.split("-")[1]

    if bond.order != "single":
        return 0.0
    if ae in {"H", "F"} or be in {"H", "F"}:
        return 0.0
    if bond.in_small_ring:
        return 0.0

    total = 0.0
    center_atoms = {bond.a, bond.b}

    for other_key, other_bond in bond_map.items():
        if other_key == bond_key:
            continue
        other_atoms = {other_bond.a, other_bond.b}
        shared = center_atoms & other_atoms
        if not shared:
            continue
        oe1 = other_bond.a.split("-")[1]
        oe2 = other_bond.b.split("-")[1]
        if other_bond.order == "single" and oe1 not in {"H", "F"} and oe2 not in {"H", "F"}:
            total += ENV_PARAMS["ESBC"]

    return total


def radical_loc_correction(atom: LocAtom) -> float:
    if atom.radical_count == 0:
        return 0.0

    total = 0.0
    for nbr in atom.neighbors:
        nbr_el = nbr.split("-")[1]
        if nbr_el == "H":
            total += ENV_PARAMS["RH"] * atom.radical_count
        else:
            total += ENV_PARAMS["RA"] * atom.radical_count

    return total


def compute_loc_correction(
    atoms: dict[str, LocAtom],
    bonds: dict[str, LocBond],
) -> LocCorrectionBreakdown:
    out = LocCorrectionBreakdown()

    for label, atom in atoms.items():
        val = atomic_loc_correction(atom)
        if abs(val) > 1e-12:
            out.atomic[label] = val

    for key, bond in bonds.items():
        val = bond_loc_correction(bond)
        if abs(val) > 1e-12:
            out.bond[key] = val

        env = environment_loc_correction(key, bond, bonds)
        if abs(env) > 1e-12:
            out.environment[key] = env

    for label, atom in atoms.items():
        val = radical_loc_correction(atom)
        if abs(val) > 1e-12:
            out.radical[label] = val

    return out
