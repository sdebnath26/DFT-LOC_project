from __future__ import annotations

import itertools
import math
from collections import defaultdict

from utils.io_utils import Molecule

# Approximate covalent radii in angstroms.
COVALENT_RADII = {
    "H": 0.31,
    "B": 0.85,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Si": 1.11,
    "Al": 1.21,
}


def distance(mol: Molecule, i: int, j: int) -> float:
    ai = mol.atoms[i]
    aj = mol.atoms[j]
    return math.dist((ai.x, ai.y, ai.z), (aj.x, aj.y, aj.z))


def angle(mol: Molecule, a: int, b: int, c: int) -> float:
    ax, ay, az = mol.atoms[a].x, mol.atoms[a].y, mol.atoms[a].z
    bx, by, bz = mol.atoms[b].x, mol.atoms[b].y, mol.atoms[b].z
    cx, cy, cz = mol.atoms[c].x, mol.atoms[c].y, mol.atoms[c].z

    p = (ax - bx, ay - by, az - bz)
    r = (cx - bx, cy - by, cz - bz)

    p_dot_r = sum(pi * ri for pi, ri in zip(p, r))
    p_norm = math.sqrt(sum(pi * pi for pi in p))
    r_norm = math.sqrt(sum(ri * ri for ri in r))

    if p_norm == 0.0 or r_norm == 0.0:
        raise ValueError("Cannot compute angle with zero-length vector")

    cos_theta = max(-1.0, min(1.0, p_dot_r / (p_norm * r_norm)))
    return math.degrees(math.acos(cos_theta))


def build_connectivity(mol: Molecule, tolerance: float = 0.40) -> dict[str, list[str]]:
    """Infer bonded neighbors by covalent radii sum + tolerance."""
    labels = [f"{i}-{atom.symbol}" for i, atom in enumerate(mol.atoms)]
    neighbors: dict[str, list[str]] = defaultdict(list)

    for i, j in itertools.combinations(range(mol.atom_count), 2):
        e1 = mol.atoms[i].symbol
        e2 = mol.atoms[j].symbol
        r1 = COVALENT_RADII.get(e1, 0.75)
        r2 = COVALENT_RADII.get(e2, 0.75)
        if distance(mol, i, j) <= (r1 + r2 + tolerance):
            neighbors[labels[i]].append(labels[j])
            neighbors[labels[j]].append(labels[i])

    # Ensure all atoms exist even if isolated.
    for label in labels:
        neighbors.setdefault(label, [])

    return dict(neighbors)


def bond_distances(mol: Molecule, connectivity: dict[str, list[str]]) -> dict[str, float]:
    out: dict[str, float] = {}
    seen: set[tuple[int, int]] = set()
    for a, nbrs in connectivity.items():
        ia = int(a.split("-")[0])
        for b in nbrs:
            ib = int(b.split("-")[0])
            key = tuple(sorted((ia, ib)))
            if key in seen:
                continue
            seen.add(key)
            out[f"{a}:{b}"] = distance(mol, ia, ib)
    return out
