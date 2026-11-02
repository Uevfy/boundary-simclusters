"""The standard run's headline numbers for one pair, as one nested dict.

`measure(sample, table)` runs every step on `sample` (the five `run()`
functions, as scripts/analyze.py does) and returns the numbers below. The
regression test (tests/golden.py, tests/golden.json) compares it with the stored
values.

The steps write their usual files to `output_paths.OUT` / `.PROCESSED`; a
caller that wants them elsewhere wraps the call in `output_paths.redirected`.
Not a pipeline step itself: it has no output of its own.
"""

from pathlib import Path

import pandas as pd

from boundary_simclusters import output_paths, samples
from boundary_simclusters.pipeline import (classify_decompositions, classify_surfaces, make_figures,
                                           match_events, paired_calopart_summary)


def counts(series):
    return {str(k): int(v) for k, v in sorted(series.items())}


def whole(value):
    """A median of counts: an int when whole ("2", not "2.00"),
    else the float (an even sample can give 32.5)."""
    value = float(value)
    return int(value) if value.is_integer() else value


def measure(sample, table):
    """Run the standard run on `sample` with SimCluster `table`; return its numbers."""
    events, caloparts = match_events.run(sample)
    paired = paired_calopart_summary.run(sample, table)
    figures = make_figures.run(sample, table)
    labels = classify_surfaces.run(sample, table)
    decomposition = classify_decompositions.run(sample, table)[1]
    unmatched = pd.read_parquet(
        output_paths.PROCESSED / f"unmatched_events_{sample.tag}.parquet")

    reordered = caloparts[caloparts.fc0_cp != caloparts.fc1_cp]
    numbers = {
        "match": {
            "events_paired": len(events),
            "caloparticles_paired": len(caloparts),
            "unmatched_fc0": int((unmatched.source_file == sample.files["fc0"].name).sum()),
            "unmatched_fc1": int((unmatched.source_file == sample.files["fc1"].name).sum()),
            "events_in_another_row": int((events.fc0_row != events.fc1_row).sum()),
            "events_caloparticles_reordered": int(reordered.fc0_row.nunique()),
        },
        "summary": {
            **{f"{stat}_{k}": value
               for k in samples.KEYS
               for stat, value in (("simclusters", int(paired[f"n_{k}"].sum())),
                                   ("mean", float(paired[f"n_{k}"].mean())),
                                   ("median", whole(paired[f"n_{k}"].median())),
                                   ("max", int(paired[f"n_{k}"].max())),
                                   ("with_zero", int((paired[f"n_{k}"] == 0).sum())),
                                   ("sum_boundary_energy", float(paired[f"sumE_{k}"].sum())))},
            "ratio_of_means": float(paired.n_fc1.mean() / paired.n_fc0.mean()),
            "share_equal_count": float((paired.n_fc0 == paired.n_fc1).mean()),
            "share_finecalo_more": float((paired.n_fc1 > paired.n_fc0).mean()),
        },
        "figures": sorted(Path(p).name for p in figures),
        "surfaces": {
            k: {
                "class": counts(labels[labels.config == k]["class"].value_counts()),
                "group": counts(labels[labels.config == k]["group"].value_counts()),
                "entering": int(labels[labels.config == k]["entering"].sum()),
            }
            for k in samples.KEYS
        },
    }
    numbers["decomposition"] = decomposition
    return numbers
