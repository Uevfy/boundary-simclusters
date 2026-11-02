"""The numbers the regression test compares, and how to regenerate them.

    python tests/golden.py            # rewrite every sample in golden.json
    python tests/golden.py e20 g20    # rewrite only these

`measure(key, workdir)` runs the standard run (the `run()` functions
`scripts/analyze.py` calls) on one local pair, writing into `workdir` instead
of outputs/ and data/processed/, and returns its numbers as a dict
(`boundary_simclusters.pipeline.numbers`).
The test (`test_standard_run.py`) and this regenerator both call it, so what is
stored and what is checked cannot drift apart.

Regenerate only when a change is *meant* to move a number, and commit the new
golden.json with that change: the diff is the record of what moved.
"""

import io
import json
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

from boundary_simclusters import output_paths, samples
from boundary_simclusters.pipeline import numbers

GOLDEN = Path(__file__).with_name("golden.json")

# tag: (Standard file, FineCalo file, table). Both files are named, so the test
# never depends on which other files sit in the folder.
SAMPLES = {
    "e20": ("bound_electron_fc0_20GeV_testNanoML.root",
            "bound_electron_fc1_20GeV_testNanoML.root", "boundary"),
    "e200": ("bound_electron_fc0_200GeV_testNanoML.root",
             "bound_electron_fc1_200GeV_testNanoML.root", "boundary"),
    "g20": ("bound_photon_fc0_20GeV_testNanoML.root",
            "bound_photon_fc1_20GeV_testNanoML.root", "boundary"),
    "g200": ("bound_photon_fc0_200GeV_testNanoML.root",
             "bound_photon_fc1_200GeV_testNanoML.root", "boundary"),
    "s1": ("testNanoML.root", "FinetestNanoML.root", "legacy"),
}


def files(key):
    """The sample's two files in the repo root."""
    fc0, fc1, _ = SAMPLES[key]
    return output_paths.REPO / fc0, output_paths.REPO / fc1


def required(key):
    """Every input of the sample that is not in git: its two files and, on the
    boundary table, the parent links its classes are built from. Without those
    the classification still runs and gives other numbers."""
    links = (output_paths.REPO / "data" / "processed" / f"ancestry_{key}_fc1.parquet",)
    return files(key) + (links if SAMPLES[key][2] == "boundary" else ())


def measure(key, workdir):
    """Run the standard run on sample `key` inside `workdir`; return its numbers."""
    with output_paths.redirected(workdir), redirect_stdout(io.StringIO()):
        sample = samples.load_sample(*files(key), tag=key)
        return numbers.measure(sample, SAMPLES[key][2])


def main(keys):
    golden = json.loads(GOLDEN.read_text()) if GOLDEN.exists() else {}
    for key in keys or SAMPLES:
        missing = [p.name for p in required(key) if not p.exists()]
        if missing:
            raise SystemExit(f"{key}: missing {', '.join(missing)}")
        with tempfile.TemporaryDirectory() as tmp:
            golden[key] = measure(key, Path(tmp))
        print(f"{key}: measured")
    GOLDEN.write_text(json.dumps({k: golden[k] for k in SAMPLES if k in golden},
                                 indent=2) + "\n")
    print(f"wrote {GOLDEN}")


if __name__ == "__main__":
    main(sys.argv[1:])
