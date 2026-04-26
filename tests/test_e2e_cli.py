"""End-to-end tests for the dft-loc CLI — run as a real subprocess.

These tests verify the user-facing behaviour of the installed CLI entry point
(or `python -m dft_loc.cli`): exit codes, stdout content, and error paths.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SRC = str(Path(__file__).resolve().parents[1] / "src")
_ENV = {**os.environ, "PYTHONPATH": _SRC}

CLI = [sys.executable, "-m", "dft_loc.cli"]

WATER_XYZ = """\
3
water
O 0.0 0.0 0.0
H 0.0 0.0 0.96
H 0.93 0.0 -0.24
"""

WATER_LOG = """\
MO # 1 Energy: -1.35 Occ: 2.0
0 O 2s   0.80
0 O 2pz  0.50
1 H 1s   0.30
MO # 2 Energy: -0.73 Occ: 2.0
0 O 2pz  0.60
1 H 1s   0.70
2 H 1s   0.70
MO # 3 Energy: -0.60 Occ: 2.0
0 O 2pz  0.40
2 H 1s   0.50
"""


def run(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [*CLI, *args],
        capture_output=True,
        text=True,
        cwd=cwd,
        env=_ENV,
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def water_files(tmp_path: Path):
    xyz = tmp_path / "water.xyz"
    log = tmp_path / "water.log"
    xyz.write_text(WATER_XYZ, encoding="utf-8")
    log.write_text(WATER_LOG, encoding="utf-8")
    return xyz, log


# ---------------------------------------------------------------------------
# E2E: no subcommand / help
# ---------------------------------------------------------------------------

class TestHelp:
    def test_no_args_exits_zero(self):
        result = run()
        assert result.returncode == 0

    def test_no_args_prints_usage(self):
        result = run()
        assert "usage" in result.stdout.lower() or "dft-loc" in result.stdout

    def test_help_flag_exits_zero(self):
        result = run("--help")
        assert result.returncode == 0

    def test_help_flag_mentions_analyze(self):
        result = run("--help")
        assert "analyze" in result.stdout

    def test_analyze_help_exits_zero(self):
        result = run("analyze", "--help")
        assert result.returncode == 0

    def test_analyze_help_shows_xyz_and_log(self):
        result = run("analyze", "--help")
        assert "--xyz" in result.stdout
        assert "--log" in result.stdout


# ---------------------------------------------------------------------------
# E2E: analyze — happy path
# ---------------------------------------------------------------------------

class TestAnalyzeHappyPath:
    def test_exits_zero(self, water_files):
        xyz, log = water_files
        result = run("analyze", "--xyz", str(xyz), "--log", str(log))
        assert result.returncode == 0

    def test_prints_bond_summary_header(self, water_files):
        xyz, log = water_files
        result = run("analyze", "--xyz", str(xyz), "--log", str(log))
        assert "Bond summary" in result.stdout

    def test_prints_valency_summary_header(self, water_files):
        xyz, log = water_files
        result = run("analyze", "--xyz", str(xyz), "--log", str(log))
        assert "Valency summary" in result.stdout

    def test_reports_atom_count(self, water_files):
        xyz, log = water_files
        result = run("analyze", "--xyz", str(xyz), "--log", str(log))
        assert "3 atoms" in result.stdout

    def test_max_occupied_mos_flag_accepted(self, water_files):
        xyz, log = water_files
        result = run(
            "analyze", "--xyz", str(xyz), "--log", str(log),
            "--max-occupied-mos", "2",
        )
        assert result.returncode == 0
        assert "Bond summary" in result.stdout

    def test_no_stderr_on_success(self, water_files):
        xyz, log = water_files
        result = run("analyze", "--xyz", str(xyz), "--log", str(log))
        assert result.stderr == ""


# ---------------------------------------------------------------------------
# E2E: analyze — error paths
# ---------------------------------------------------------------------------

class TestAnalyzeErrors:
    def test_missing_xyz_arg_exits_nonzero(self, water_files):
        _, log = water_files
        result = run("analyze", "--log", str(log))
        assert result.returncode != 0

    def test_missing_log_arg_exits_nonzero(self, water_files):
        xyz, _ = water_files
        result = run("analyze", "--xyz", str(xyz))
        assert result.returncode != 0

    def test_nonexistent_xyz_file_raises(self, water_files, tmp_path):
        _, log = water_files
        result = run("analyze", "--xyz", str(tmp_path / "no.xyz"), "--log", str(log))
        assert result.returncode != 0

    def test_nonexistent_log_file_raises(self, water_files, tmp_path):
        xyz, _ = water_files
        result = run("analyze", "--xyz", str(xyz), "--log", str(tmp_path / "no.log"))
        assert result.returncode != 0

    def test_malformed_xyz_exits_nonzero(self, tmp_path, water_files):
        _, log = water_files
        bad_xyz = tmp_path / "bad.xyz"
        bad_xyz.write_text("not_a_number\ntitle\nH 0 0 0\n", encoding="utf-8")
        result = run("analyze", "--xyz", str(bad_xyz), "--log", str(log))
        assert result.returncode != 0

    def test_unknown_subcommand_exits_nonzero(self):
        result = run("no_such_command")
        assert result.returncode != 0
