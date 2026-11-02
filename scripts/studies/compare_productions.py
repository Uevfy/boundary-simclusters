"""One gun sample from two productions: does anything change between them?

    python scripts/studies/compare_productions.py \\
        bound_electron_fc0_20GeV_testNanoML.root bound_electron_fc1_20GeV_testNanoML.root \\
        d121/rerun_pre1eta_electron_fc0_20GeV_testNanoML.root \\
        d121/rerun_pre1eta_electron_fc1_20GeV_testNanoML.root \\
        --root d121 --tag e20

Old Standard, old FineCalo, new Standard, new FineCalo, in that order. Built for
master vs pre1_D121 (`d121/`). What it compares, nearest the file first:

    overlap    events both productions hold, per config, and where the others sit
    branches   every branch both files carry, event by event on the GEN key.
               A branch that differs is compared again with each event's values
               sorted, which separates a reordering from a change of value
    outcome    the standard run's result per CaloParticle (counts),
               on the CaloParticles of events all four files hold
    headline   `pipeline.numbers.measure` on each pair, side by side

Files are joined on the GEN key and CaloParticles on (pdgId, eta, phi), never by
index or track ID. Both standard runs go into a temporary folder, so an old
parquet cannot enter.

Writes outputs/comparison/report_production_<tag>.txt and
data/processed/production_comparison_<tag>.parquet (one row per branch and config).
"""

import argparse
import io
import sys
import tempfile
from collections import Counter
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas as pd

from boundary_simclusters import output_paths, samples
from boundary_simclusters.matching import genpart_signature
from boundary_simclusters.pipeline import numbers

OUTCOME_COLUMNS = ["n_std", "n_fine", "n_shared", "n_side_entries"]  # class needs each production's ancestry


def same_values(a, b):
    """Elementwise equality per event of two jagged columns, NaN equal to NaN."""
    if len(a) == 0:
        return np.ones(0, dtype=bool)
    if a.dtype != object:                  # a flat branch: one value per event
        same = a == b
        if a.dtype.kind == "f":
            same |= np.isnan(a) & np.isnan(b)
        return same
    out = np.zeros(len(a), dtype=bool)
    for i, (x, y) in enumerate(zip(a, b)):
        if len(x) != len(y):
            continue
        equal = x == y
        if x.dtype.kind == "f":
            equal = equal | (np.isnan(x) & np.isnan(y))
        out[i] = bool(np.all(equal))
    return out


def same_when_sorted(a, b):
    """The same, ignoring the order of the values inside each event."""
    if a.dtype != object:
        return same_values(a, b)
    def sorted_events(column):
        out = np.empty(len(column), dtype=object)
        out[:] = [np.sort(x) for x in column]
        return out
    return same_values(sorted_events(a), sorted_events(b))


def event_rows(old, new):
    """Rows of the events both trees hold, on GEN keys occurring once in each."""
    keys_old, keys_new = genpart_signature(old), genpart_signature(new)
    count_old, count_new = Counter(keys_old), Counter(keys_new)
    shared = {k for k in count_old if count_old[k] == 1 and count_new[k] == 1}
    row_new = {k: r for r, k in enumerate(keys_new) if k in shared}
    pairs = [(r, row_new[k]) for r, k in enumerate(keys_old) if k in shared]
    only_old = [r for r, k in enumerate(keys_old) if k not in shared]
    only_new = [r for r, k in enumerate(keys_new) if k not in shared]
    rows_old, rows_new = (np.array(x, dtype=int) for x in zip(*pairs))
    return rows_old, rows_new, only_old, only_new


def compare_overlap(key, old, new):
    rows_old, rows_new, only_old, only_new = event_rows(old, new)
    print(f"  {samples.LABELS[key]:9s} shared {len(rows_old)} of {old.num_entries} old / "
          f"{new.num_entries} new; in the same row: {int((rows_old == rows_new).sum())}")
    print(f"            only old: {len(only_old)} at rows {only_old}")
    print(f"            only new: {len(only_new)} at rows {only_new}")
    return rows_old, rows_new


def compare_branches(key, old, new, rows_old, rows_new):
    """One row per branch both trees carry: events that differ, as stored and sorted."""
    records = []
    for branch in sorted(set(old.keys()) & set(new.keys())):
        a = old[branch].array(library="np")[rows_old]
        b = new[branch].array(library="np")[rows_new]
        same = same_values(a, b)
        differing = int((~same).sum())
        unsorted_only = 0
        if differing:
            unsorted_only = differing - int((~same_when_sorted(a[~same], b[~same])).sum())
        records.append({"config": key, "branch": branch, "events": len(a),
                        "events_differing": differing, "reordered_only": unsorted_only})
    return pd.DataFrame(records)


def print_branches(table):
    for key, group in table.groupby("config", sort=False):
        changed = group[group.events_differing > 0]
        print(f"  {samples.LABELS[key]:9s} {len(group)} branches, "
              f"{len(group) - len(changed)} identical on every shared event, {len(changed)} differ")
        if changed.empty:
            continue
        collections = changed.branch.str.split("_").str[0]
        for name, part in changed.groupby(collections):
            real = part[part.events_differing > part.reordered_only]
            print(f"    {name}: {len(part)} branches, up to {part.events_differing.max()} of "
                  f"{part.events.iloc[0]} events; order only: {len(part) - len(real)}, "
                  f"values change: {len(real)}")
            for r in real.sort_values("events_differing", ascending=False).head(5).itertuples():
                print(f"      {r.branch}: {r.events_differing} events")


def run_standard(sample):
    """`numbers.measure` in a temporary folder, plus the per-CaloParticle tables."""
    with tempfile.TemporaryDirectory() as tmp, output_paths.redirected(tmp), \
            redirect_stdout(io.StringIO()):
        measured = numbers.measure(sample, sample.table)
        processed = output_paths.PROCESSED
        classes = pd.read_parquet(processed / f"decompositions_{sample.tag}_bound.parquet")
        paired = pd.read_parquet(processed / f"paired_calopart_{sample.tag}_bound.parquet")
    return measured, classes, paired


def calopart_outcome(sample, classes, paired):
    """Per CaloParticle: (GEN key of its event, pdgId, eta, phi) -> outcome."""
    gen = genpart_signature(sample.tree("fc0"))
    cps = samples.caloparticle_frame(sample.tree("fc0")).set_index(["event", "cp"])
    rows = classes.set_index(["fc0_row", "fc0_cp"])
    paired = paired.set_index(["fc0_row", "calopart"])
    records = {}
    for (row, cp), r in rows.iterrows():
        c = cps.loc[(row, cp)]
        p = paired.loc[(row, cp)]
        key = (gen[row], int(c.pdg), float(c.eta), float(c.phi))
        records[key] = {**{col: r[col] for col in OUTCOME_COLUMNS},
                        "n_fc0": int(p.n_fc0), "n_fc1": int(p.n_fc1),
                        "sumE_fc0": float(p.sumE_fc0), "sumE_fc1": float(p.sumE_fc1)}
    return records


def compare_outcome(old, new):
    shared = old.keys() & new.keys()
    print(f"  CaloParticles: old {len(old)}, new {len(new)}, in both {len(shared)}")
    fields = list(next(iter(old.values())))
    for field in fields:
        n = sum(old[k][field] != new[k][field] for k in shared)
        print(f"    {field:12s} differs for {n}")
    changed = [k for k in shared if old[k] != new[k]]
    print(f"  CaloParticles with any difference: {len(changed)}")
    for k in changed[:10]:
        print(f"    pdg {k[1]} eta {k[2]:.3f}: "
              + ", ".join(f"{f} {old[k][f]} -> {new[k][f]}" for f in fields if old[k][f] != new[k][f]))


def leaves(tree, path=""):
    if isinstance(tree, dict):
        for k, v in tree.items():
            yield from leaves(v, f"{path}.{k}" if path else k)
    else:
        yield path, tree


def print_headline(old, new):
    a, b = dict(leaves(old)), dict(leaves(new))
    moved = sorted(k for k in a.keys() | b.keys() if a.get(k) != b.get(k) and k != "figures")
    print(f"  {len(a)} numbers, {len(moved)} move")
    for k in moved:
        print(f"    {k:42s} {a.get(k)}  ->  {b.get(k)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs=4, type=Path, metavar="ROOT",
                    help="old fc0, old fc1, new fc0, new fc1")
    ap.add_argument("--root", type=Path, required=True, help="tree the results go to (d121)")
    ap.add_argument("--tag", required=True, help="sample label, e.g. e20")
    args = ap.parse_args()

    old = samples.load_sample(args.files[0], args.files[1], tag=args.tag)
    new = samples.load_sample(args.files[2], args.files[3], tag=f"{args.tag}_new")
    if (old.particle, old.energy) != (new.particle, new.energy):
        raise SystemExit(f"not one gun sample: {old.particle} {old.energy} vs {new.particle} {new.energy}")

    with output_paths.rooted(args.root):
        report_dir = output_paths.OUT / "comparison"
        report_dir.mkdir(parents=True, exist_ok=True)
        output_paths.PROCESSED.mkdir(parents=True, exist_ok=True)
        table_path = output_paths.PROCESSED / f"production_comparison_{args.tag}.parquet"
    report = io.StringIO()
    with redirect_stdout(report):
        print(f"{old.particle} {old.energy}")
        for label, s in (("old", old), ("new", new)):
            print(f"  {label}: {s.files['fc0'].name}, {s.files['fc1'].name}")

        added = sorted(set(new.tree("fc0").keys()) - set(old.tree("fc0").keys()))
        removed = sorted(set(old.tree("fc0").keys()) - set(new.tree("fc0").keys()))
        print(f"\nbranches: {len(removed)} removed, {len(added)} added (not compared)")
        for name, n in pd.Series([b.split("_")[0] for b in added]).value_counts().sort_index().items():
            print(f"  {name}: {n}")

        print("\noverlap")
        rows = {k: compare_overlap(k, old.tree(k), new.tree(k)) for k in samples.KEYS}

        print("\nbranches on the shared events")
        table = pd.concat([compare_branches(k, old.tree(k), new.tree(k), *rows[k])
                           for k in samples.KEYS], ignore_index=True)
        table.to_parquet(table_path, index=False)
        print_branches(table)

        print("\nstandard run per CaloParticle")
        numbers_old, classes_old, paired_old = run_standard(old)
        numbers_new, classes_new, paired_new = run_standard(new)
        compare_outcome(calopart_outcome(old, classes_old, paired_old),
                        calopart_outcome(new, classes_new, paired_new))

        print("\nheadline numbers, old -> new")
        print_headline(numbers_old, numbers_new)

    text = report.getvalue()
    (report_dir / f"report_production_{args.tag}.txt").write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
