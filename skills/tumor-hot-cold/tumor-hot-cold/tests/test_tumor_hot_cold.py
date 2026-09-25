"""
Tests for the Tumor Hot Cold skill.

Run:
    pytest skills/tumor-hot-cold/tests/ -v
or via ClawBio runner:
    python -m pytest skills/tumor-hot-cold/tests/ -v
"""

import json
from pathlib import Path

import pytest
import sys

# Make the skill importable regardless of working directory
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tumor_hot_cold import run, run_demo


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_output(tmp_path: Path) -> Path:
    """Provide a temporary output directory per test."""
    return tmp_path / "output"


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_demo_runs(tmp_output: Path) -> None:
    """Demo mode should complete without raising an exception."""
    run_demo(tmp_output)


def test_report_generated(tmp_output: Path) -> None:
    """Demo should produce a non-empty report.md."""
    run_demo(tmp_output)
    report = tmp_output / "report.md"
    assert report.exists(), "report.md was not created"
    assert report.stat().st_size > 0, "report.md is empty"


def test_result_json_valid(tmp_output: Path) -> None:
    """result.json should be present and parse as valid JSON."""
    run_demo(tmp_output)
    result_file = tmp_output / "result.json"
    assert result_file.exists(), "result.json was not created"
    data = json.loads(result_file.read_text())
    assert isinstance(data, dict), "result.json should be a JSON object"
    assert data.get("skill") == "tumor-hot-cold", "result.json should include skill name"


def test_reproducibility_bundle(tmp_output: Path) -> None:
    """Reproducibility bundle must include commands, environment, and checksums."""
    run_demo(tmp_output)
    repro = tmp_output / "reproducibility"
    assert (repro / "commands.sh").exists(), "commands.sh missing"
    assert (repro / "environment.yml").exists(), "environment.yml missing"
    assert (repro / "checksums.sha256").exists(), "checksums.sha256 missing"


# ---------------------------------------------------------------------------
# TODO: Add domain-specific tests below
# ---------------------------------------------------------------------------

# def test_core_logic(tmp_output: Path) -> None:
#     """TODO: Test the actual analysis logic with known input/output."""
#     demo_input = tmp_output / "test_input.txt"
#     demo_input.parent.mkdir(parents=True, exist_ok=True)
#     demo_input.write_text("...")
#     results = run(demo_input, tmp_output)
#     assert results["..."] == expected_value
