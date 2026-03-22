from __future__ import annotations

import re
from collections import defaultdict

MO_HEADER = re.compile(r"MO\s+#?\s*(\d+).*(Occ:\s*([0-9.]+))")
FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def parse_mo_coefficients(lines: list[str], occupied_threshold: float = 1.5) -> dict[int, dict[str, float]]:
    """
    Parse occupied MO coefficient sections from quantum-chemistry logs.

    Expected shape (loosely):
      MO # 3 Energy: ... Occ: 2
      0 C 2pz    0.1234
      1 O 2pz    0.5678
    """
    mo_dict: dict[int, dict[str, float]] = defaultdict(dict)
    current_mo: int | None = None
    collect = False

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        header_match = MO_HEADER.search(line)
        if header_match:
            mo_idx = int(header_match.group(1))
            occ = float(header_match.group(3))
            collect = occ >= occupied_threshold
            current_mo = mo_idx if collect else None
            continue

        if current_mo is None:
            continue

        parts = line.split()
        if len(parts) < 4:
            continue

        # require leading integer atom index
        if not parts[0].isdigit():
            continue

        atom_index = parts[0]
        symbol = parts[1]
        coeff_match = FLOAT_RE.search(parts[-1])
        if coeff_match is None:
            continue

        atom_label = f"{atom_index}-{symbol}"
        coeff = float(coeff_match.group(0))
        mo_dict[current_mo][atom_label] = mo_dict[current_mo].get(atom_label, 0.0) + abs(coeff)

    return dict(mo_dict)


def build_atom_mo_coefficients(
    mo_dict: dict[int, dict[str, float]],
    max_occupied_mos: int | None = None,
) -> dict[str, list[tuple[int, float]]]:
    """Convert MO-centric coefficients to atom-centric sorted lists."""
    atom_map: dict[str, list[tuple[int, float]]] = defaultdict(list)

    mo_levels = sorted(mo_dict)
    if max_occupied_mos is not None:
        mo_levels = mo_levels[-max_occupied_mos:]

    for mo in mo_levels:
        for atom_label, coeff in mo_dict[mo].items():
            atom_map[atom_label].append((mo, coeff))

    for atom_label in atom_map:
        atom_map[atom_label].sort(key=lambda item: item[1], reverse=True)

    return dict(atom_map)


def get_top_atom_coeff(
    atom_coeff: dict[str, list[tuple[int, float]]],
    central_atom: str | None = None,
    num_coeff: int = 4,
) -> dict[str, list[tuple[int, float]]]:
    out: dict[str, list[tuple[int, float]]] = {}
    for atom, pairs in atom_coeff.items():
        if atom == central_atom:
            out[atom] = pairs[: max(1, num_coeff)]
        elif atom.endswith("-H"):
            out[atom] = pairs[:1]
        elif atom.endswith("-B") or atom.endswith("-Al"):
            out[atom] = pairs[:3]
        else:
            out[atom] = pairs[:4]
    return out
