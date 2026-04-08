from __future__ import annotations

from collections import defaultdict

from utils.io_utils import element_symbols, atom_labels, Molecule
from analysis.loc_corrections import LocAtom, LocBond


def build_neighbor_map(bond_types: dict[str, str]) -> dict[str, list[str]]:
    nbrs = defaultdict(list)
    for pair, order in bond_types.items():
        a, b = pair.split(":")
        nbrs[a].append(b)
        nbrs[b].append(a)
    return dict(nbrs)


def infer_hybridization(
    label: str,
    element: str,
    bond_types: dict[str, str],
    neighbor_map: dict[str, list[str]],
) -> str | None:
    orders = []
    for pair, order in bond_types.items():
        a, b = pair.split(":")
        if label == a or label == b:
            orders.append(order)

    n_double = sum(1 for x in orders if x == "double")
    n_triple = sum(1 for x in orders if x == "triple")
    coordination = len(orders)

    if n_triple >= 1:
        return "sp"
    if n_double >= 1:
        return "sp2"
    if coordination >= 4:
        return "sp3"
    if coordination == 3:
        return "sp2"
    if coordination <= 2:
        return "sp"
    return None


def assign_loc_atoms(mol: Molecule, bond_types: dict[str, str]) -> dict[str, LocAtom]:
    labels = atom_labels(mol)
    elems = element_symbols(mol)
    neighbor_map = build_neighbor_map(bond_types)

    out = {}
    for label, element in zip(labels, elems):
        hyb = infer_hybridization(label, element, bond_types, neighbor_map)
        out[label] = LocAtom(
            label=label,
            element=element,
            neighbors=tuple(sorted(neighbor_map.get(label, []))),
            hybridization=hyb,
        )
    return out


def assign_loc_bonds(
    bond_types: dict[str, str],
    distances: dict[str, float],
) -> dict[str, LocBond]:
    out = {}

    for pair, order in bond_types.items():
        a, b = pair.split(":")
        dist = distances.get(pair, None)
        ae = a.split("-")[1]
        be = b.split("-")[1]

        bond = LocBond(a=a, b=b, order=order)

        # H-containing bonds
        if "H" in (ae, be):
            other = be if ae == "H" else ae
            if other in {"O", "F", "Cl"}:
                bond = LocBond(a=a, b=b, order=order, polarized_class="POLH")
            elif other == "S":
                bond = LocBond(a=a, b=b, order=order, polarized_class="IPOLH")
            else:
                bond = LocBond(a=a, b=b, order=order, polarized_class="NPOLH")

        # F-containing bonds
        elif "F" in (ae, be):
            other = be if ae == "F" else ae
            if order == "single":
                if other == "N":
                    bond = LocBond(a=a, b=b, order=order, polarized_class="IPOLF")
                else:
                    bond = LocBond(a=a, b=b, order=order, polarized_class="POLF")

        # Non-H/F bonds
        else:
            if order == "double":
                bond = LocBond(a=a, b=b, order=order, multiple_class="DBC")
            elif order == "triple":
                if ae == be:
                    bond = LocBond(a=a, b=b, order=order, multiple_class="TBNPOL")
                else:
                    bond = LocBond(a=a, b=b, order=order, multiple_class="TBPOL")
            elif order == "single":
                # very simple first pass; refine later with ring and CT logic
                if dist is not None and dist < 1.45:
                    cls = "SSBC"
                elif dist is not None and dist < 1.85:
                    cls = "MSBC"
                else:
                    cls = "LSBC"
                bond = LocBond(a=a, b=b, order=order, length_class=cls)

        out[pair] = bond

    return out
