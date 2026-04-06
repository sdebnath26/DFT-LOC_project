"""Basic import and smoke tests for dft-loc."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest


EXAMPLES_DIR = Path(__file__).parent.parent / "examples"


def test_package_imports():
    """Ensure the top-level package and all sub-modules import without error."""
    import dft_loc  # noqa: F401
    from dft_loc.utils import io_utils, geometry_utils, mo_utils, print_utils  # noqa: F401
    from dft_loc.analysis import bond_analysis  # noqa: F401
    from dft_loc import cli  # noqa: F401


def test_version_attribute():
    import dft_loc
    assert hasattr(dft_loc, "__version__")
    assert isinstance(dft_loc.__version__, str)


def test_open_xyz():
    from dft_loc.utils.io_utils import open_xyz

    xyz_path = EXAMPLES_DIR / "pyrrole.xyz"
    if not xyz_path.exists():
        pytest.skip("Example XYZ file not present")

    mol = open_xyz(xyz_path)
    assert mol.atom_count > 0
    assert len(mol.atoms) == mol.atom_count


def test_open_log():
    from dft_loc.utils.io_utils import open_log

    log_path = EXAMPLES_DIR / "fb_pyrrole.log"
    if not log_path.exists():
        pytest.skip("Example log file not present")

    lines = open_log(log_path)
    assert isinstance(lines, list)
    assert len(lines) > 0


def test_full_pipeline():
    """Run the entire analysis pipeline on the example files."""
    from dft_loc.utils.io_utils import open_xyz, open_log, atom_labels
    from dft_loc.utils.geometry_utils import bond_distances, build_connectivity
    from dft_loc.utils.mo_utils import (
        build_atom_mo_coefficients,
        get_top_atom_coeff,
        parse_mo_coefficients,
    )
    from dft_loc.analysis.bond_analysis import (
        classify_bonds,
        infer_bond_pair_indices,
        valency_atom,
    )

    xyz_path = EXAMPLES_DIR / "pyrrole.xyz"
    log_path = EXAMPLES_DIR / "fb_pyrrole.log"
    if not xyz_path.exists() or not log_path.exists():
        pytest.skip("Example files not present")

    mol = open_xyz(xyz_path)
    log_lines = open_log(log_path)

    mo_dict = parse_mo_coefficients(log_lines)
    atom_coeff = build_atom_mo_coefficients(mo_dict)
    top_coeff = get_top_atom_coeff(atom_coeff)

    connectivity = build_connectivity(mol)
    distances = bond_distances(mol, connectivity)
    bond_indices = infer_bond_pair_indices(top_coeff, connectivity)
    bond_types = classify_bonds(bond_indices)

    assert len(bond_types) > 0
    for pair, btype in bond_types.items():
        assert btype in {"single", "double", "triple", "unknown"}

    for label in atom_labels(mol):
        v = valency_atom(label, bond_types)
        assert v >= 0


def test_cli_help(capsys):
    """CLI --help should print usage and exit 0."""
    from dft_loc.cli import main

    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "dft-loc" in captured.out


def test_cli_analyze(capsys):
    """CLI analyze subcommand should produce output."""
    from dft_loc.cli import main

    xyz_path = EXAMPLES_DIR / "pyrrole.xyz"
    log_path = EXAMPLES_DIR / "fb_pyrrole.log"
    if not xyz_path.exists() or not log_path.exists():
        pytest.skip("Example files not present")

    main(["analyze", "--xyz", str(xyz_path), "--log", str(log_path)])
    captured = capsys.readouterr()
    assert "Bond summary" in captured.out
    assert "Valency summary" in captured.out
