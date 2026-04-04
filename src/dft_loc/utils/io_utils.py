from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class AtomRecord:
    symbol: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class Molecule:
    atom_count: int
    title: str
    atoms: tuple[AtomRecord, ...]


@dataclass(frozen=True)
class OrbitalContribution:
    atom_label: str
    coefficient: float


def open_xyz(path: str | Path, encoding: str = "utf-8") -> Molecule:
    """Read an XYZ file into a structured Molecule object."""
    text = Path(path).read_text(encoding=encoding).strip()
    rows = [row for row in text.splitlines() if row.strip()]
    if len(rows) < 2:
        raise ValueError(f"Invalid XYZ file (too short): {path}")

    atom_count = int(rows[0].strip())
    title = rows[1].strip()
    atom_rows = rows[2:]

    # Some files include an extra charge/multiplicity line after the title.
    if len(atom_rows) == atom_count + 1:
        atom_rows = atom_rows[1:]

    if len(atom_rows) != atom_count:
        raise ValueError(
            f"Invalid XYZ file {path}: expected {atom_count} atom lines, found {len(atom_rows)}"
        )

    atoms: list[AtomRecord] = []
    for row in atom_rows:
        parts = row.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ atom row: {row}")
        atoms.append(AtomRecord(parts[0], float(parts[1]), float(parts[2]), float(parts[3])))

    return Molecule(atom_count=atom_count, title=title, atoms=tuple(atoms))


def open_log(path: str | Path, encoding: str = "utf-8") -> list[str]:
    """Return raw log lines."""
    return Path(path).read_text(encoding=encoding).splitlines()


def atom_labels(mol: Molecule) -> tuple[str, ...]:
    """Return labels in index-symbol format used by legacy scripts (e.g. '0-C')."""
    return tuple(f"{i}-{atom.symbol}" for i, atom in enumerate(mol.atoms))


def element_symbols(mol: Molecule) -> tuple[str, ...]:
    return tuple(atom.symbol for atom in mol.atoms)
