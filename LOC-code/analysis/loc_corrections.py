from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

BondType = Literal["single", "double", "triple", "unknown"]

# ----------------------------
# Parameter tables from LOC paper
# B3LYP-LOC values fit to full G3 set
# ----------------------------

ATOM_PARAMS = {
    "Be_sp": 7.31,
    "N_sp": 0.00,
    "N_sp2": 4.46,
    "N_sp3": 3.68,
    "N_quart": 6.04,
    "P_sp": 0.00,
    "P_sp2": 4.46,
    "P_sp3": 3.68,
    "P_quart": 6.04,
    "O_sp2": 0.95,
    "O_sp3": 1.78,
    "OCT_EXP": 4.84,
}

BOND_PARAMS = {
    "NPOLH": 0.37,
    "IPOLH": 0.00,   # S-H
    "POLH": -1.42,

    "NPOLF": 0.88,
    "IPOLF": 0.00,   # N-F
    "POLF": -0.95,

    "SSBC": -1.41,
    "MSBC": -1.96,
    "LSBC": -2.63,

    "DBC": -1.03,
    "TBNPOL": -1.68,
    "TBPOL": 0.88,

    "CT": -4.53,
}

ENV_PARAMS = {
    "ESBC": -0.49,
    "RH": 0.52,
    "RA": 1.65,
    "RT": -2.58,
}


@dataclass(frozen=True)
class LocAtom:
    label: str
    element: str
    neighbors: tuple[str, ...]
    formal_charge: int = 0
    radical_count: int = 0
    hybridization: str | None = None
    octet_expanded: bool = False


@dataclass(frozen=True)
class LocBond:
    a: str
    b: str
    order: BondType
    aromatic_like: bool = False
    in_small_ring: bool = False
    charge_transfer: bool = False
    polarized_class: str | None = None   # NPOLH, POLH, ...
    length_class: str | None = None      # SSBC, MSBC, LSBC
    multiple_class: str | None = None    # DBC, TBNPOL, TBPOL


@dataclass
class LocCorrectionBreakdown:
    atomic: dict[str, float] = field(default_factory=dict)
    bond: dict[str, float] = field(default_factory=dict)
    environment: dict[str, float] = field(default_factory=dict)
    radical: dict[str, float] = field(default_factory=dict)

    @property
    def total(self) -> float:
        return (
            sum(self.atomic.values())
            + sum(self.bond.values())
            + sum(self.environment.values())
            + sum(self.radical.values())
        )
