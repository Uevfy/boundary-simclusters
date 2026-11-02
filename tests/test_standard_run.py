"""Regression test: the standard run still gives the numbers in golden.json.

    pytest                                            # every sample, ~1 min
    pytest -k g200                                    # samples whose name contains g200
    pytest "tests/test_standard_run.py::test_standard_run[e20]"   # exactly e20

(`-k` matches substrings, so `-k e20` runs e200 as well.)

Each sample runs the pipeline steps into a temporary folder (outputs/ and
data/processed/ are never touched) and compares every stored number: counts
exactly, floats to a relative 1e-6. A sample whose inputs are not there (ROOT
files in the repo root, ancestry parquet in data/processed/) is skipped. What is measured, and how to regenerate golden.json when
a change is meant to move a number, is in golden.py.
"""

import json
import math

import pytest

import golden

STORED = json.loads(golden.GOLDEN.read_text())
REL_TOL = 1e-6


def differences(expected, got, path=""):
    """Every leaf where `got` departs from `expected`, as readable lines."""
    if isinstance(expected, dict) and isinstance(got, dict):
        lines = []
        for key in sorted(set(expected) | set(got)):
            where = f"{path}.{key}" if path else key
            if key not in got:
                lines.append(f"{where}: expected {expected[key]!r}, missing")
            elif key not in expected:
                lines.append(f"{where}: not in golden.json, got {got[key]!r}")
            else:
                lines += differences(expected[key], got[key], where)
        return lines
    if isinstance(expected, float) or isinstance(got, float):
        same = math.isclose(expected, got, rel_tol=REL_TOL, abs_tol=1e-12)
    else:
        same = expected == got
    return [] if same else [f"{path}: expected {expected!r}, got {got!r}"]


@pytest.mark.parametrize("key", list(golden.SAMPLES))
def test_standard_run(key, tmp_path):
    missing = [p.name for p in golden.required(key) if not p.exists()]
    if missing:
        pytest.skip(f"input files not found: {', '.join(missing)}")
    if key not in STORED:
        pytest.fail(f"no golden numbers for {key}; run: python tests/golden.py {key}")
    lines = differences(STORED[key], golden.measure(key, tmp_path))
    assert not lines, f"{key} moved ({len(lines)}):\n  " + "\n  ".join(lines)
