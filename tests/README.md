# tests

A golden-number regression test for the standard run. It reruns the steps
of `boundary_simclusters.pipeline` (the ones `scripts/analyze.py` calls) on
five local pairs and checks that every number still matches `golden.json`.

```bash
pytest                                                        # all five, ~1 min
pytest "tests/test_standard_run.py::test_standard_run[e20]"   # one sample
```

| sample | table | why it is in the test |
|---|---|---|
| `e20` | boundary | the quickest full run |
| `e200` | boundary | electrons at 200 GeV |
| `g20` | boundary | photons; 2 events whose CaloParticles are stored in another order (docs/data_handling.md) |
| `g200` | boundary | photons at 200 GeV; 1 reordered event |
| `s1` | legacy | the legacy collection (995 events, 5.35×), and the only `interior` SimClusters among the five |

**What is compared**, per sample (`golden.py`, `measure`):
- **match:** events and CaloParticles paired, unmatched rows, events in another row, events with reordered CaloParticles
- **summary:** SimClusters per config (total, mean, max, CaloParticles with none), Σ boundaryEnergy, ratio of means, share with equal counts, share where FineCalo has more
- **figures:** the names of the files written. They are not compared by pixel, so a deliberate change of style passes.
- **surfaces:** SimClusters per config by class and by face group, and how many are entering
- **decomposition** (boundary samples only): CaloParticles per class, the track-ID
  agreement, shared, replaced and descending SimClusters.

Counts must match exactly, floats to a relative 1e-6. A failure lists every
number that moved, e.g. `summary.ratio_of_means: expected 2.38…, got 2.41…`.

**Where it writes.** Each sample runs in a temporary folder: `output_paths.OUT`
and `output_paths.PROCESSED` are pointed there for the run, so `outputs/` and
`data/processed/` are never touched. Library code must therefore read those two
as `output_paths.OUT` / `output_paths.PROCESSED` when it writes, never copy them
at import.

**Missing data.** The inputs are not in git: the ROOT files in the repo root and, for
the boundary samples, `data/processed/ancestry_<tag>_fc1.parquet` (written from RECO by
`scripts/tools/extract_ancestry.py`). A sample missing any of them is skipped, and the
skip says which file is missing. `golden.py` refuses to regenerate such a sample.

**When a number is meant to move**, for example after fixing a bug or changing a
definition:

```bash
python tests/golden.py          # all samples
python tests/golden.py g20      # only this one
```

Commit the new `golden.json` together with the change that moved it, and say
why in the commit message: the diff is the record of what moved. Never
regenerate just to make a failing test pass without knowing why it failed.
The test proves that the code still gives what it gave, not that it is right.
