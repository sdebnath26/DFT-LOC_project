from __future__ import annotations

import argparse
from pathlib import Path

from analysis.bond_analysis import (
    bond_order1,
    classify_bonds,
    infer_bond_pair_indices,
    valency_atom,
)
from utils.geometry_utils import bond_distances, build_connectivity
from utils.io_utils import atom_labels, open_log, open_xyz
from utils.mo_utils import build_atom_mo_coefficients, get_top_atom_coeff, parse_mo_coefficients


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Localized-orbital bond analysis from XYZ geometry + MO log output"
    )
    parser.add_argument("--xyz", default="data/pyrrole.xyz", help="Path to XYZ file")
    parser.add_argument("--log", default="data/fb_pyrrole.log", help="Path to log file")
    parser.add_argument(
        "--max-occupied-mos",
        type=int,
        default=None,
        help="Use only the last N occupied MOs for bond inference",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    xyz_path = Path(args.xyz)
    log_path = Path(args.log)

    mol = open_xyz(xyz_path)
    log_lines = open_log(log_path)

    mo_dict = parse_mo_coefficients(log_lines)
    atom_coeff = build_atom_mo_coefficients(mo_dict, max_occupied_mos=args.max_occupied_mos)
    top_coeff = get_top_atom_coeff(atom_coeff)

    connectivity = build_connectivity(mol)
    distances = bond_distances(mol, connectivity)
    bond_indices = infer_bond_pair_indices(top_coeff, connectivity)
    bond_types = classify_bonds(bond_indices)

    print(f"Loaded: {xyz_path} ({mol.atom_count} atoms)")
    print(f"Loaded: {log_path} ({len(log_lines)} lines)")
    print("\nBond summary")
    for pair in sorted(bond_types):
        dist = distances.get(pair)
        dist_text = f"{dist:.3f} Å" if dist is not None else "n/a"
        print(f"- {pair:20s} {bond_types[pair]:7s}  shared_mos={bond_indices[pair]}  d={dist_text}")

    print("\nValency summary")
    for label in atom_labels(mol):
        v = valency_atom(label, bond_types)
        s, d, t = bond_order1(label, bond_types)
        print(f"- {label:6s} valency={v} (single={s}, double={d}, triple={t})")


if __name__ == "__main__":
    main()
