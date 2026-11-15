"""The standard run on one fc0/fc1 pair, in one command: nanoML file(s) in, results out.

    python scripts/analyze.py bound_electron_fc0_20GeV_testNanoML.root
    python scripts/analyze.py a.root b.root --table legacy --tag s1
    python scripts/analyze.py a.root --steps summary figures
    python scripts/analyze.py a.root --steps figures --only energy n_per_event
    python scripts/analyze.py other/a.root other/b.root --root other --tag e20_other

Give one file and its partner is found in the same folder, or give both in
either order; which is Standard and which FineCalo, the particle, the energy
and the tables are read from the contents (`samples.py`). The steps of
`boundary_simclusters.pipeline` then run in order:

    match      event and CaloParticle mappings       data/processed/{event,calopart}_matching_<tag>.parquet
    summary    paired per-CaloParticle comparison   data/processed/paired_calopart_<tag>[_bound].parquet
    figures    the standard figures                 outputs/<particle>/<energy>/<table>/
    surfaces   which boundary face each crossing hit data/processed/surface_<tag>[_bound].parquet + figures
    decompositions  each CaloParticle's decomposition    data/processed/decompositions_<tag>_bound.parquet
               compared by root track (each table its own classes;
               parent links from scripts/tools/extract_ancestry.py)
    compare    the figures and surface figures with  outputs/<particle>/<energy>/comparisons/
               both collections overlaid             (pairs carrying both tables only)

Everything printed is also written to report_<tag>[_bound].txt in the figure
folder, and the run ends with the headline numbers. Boundary table by default
where both files carry it, `--table legacy` for the legacy one. `--steps` runs a subset (summary needs the mappings of a match run),
`--only` a subset of the figures.
"""

import argparse
import sys
from contextlib import redirect_stdout
from pathlib import Path

from boundary_simclusters import output_paths, samples
from boundary_simclusters.output_paths import table_dir
from boundary_simclusters.pipeline import (classify_decompositions, classify_surfaces, make_figures,
                                           match_events, paired_calopart_summary)

STEPS = ("match", "summary", "figures", "surfaces", "decompositions", "compare")


class Tee:
    """Write to several streams at once: the terminal and the report file."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for s in self.streams:
            s.write(text)

    def flush(self):
        for s in self.streams:
            s.flush()


def banner(title):
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


def headline(sample, table, results):
    """The few numbers worth reading first; the full report is above."""
    banner("Summary")
    print(f"sample   {sample.tag}: {sample.particle} {sample.energy}, {table} table")
    for key in samples.KEYS:
        print(f"  {samples.LABELS[key]:9s} {sample.files[key].name}")
    if "match" in results:
        events, caloparts = results["match"]
        print(f"paired   {len(events)} events, {len(caloparts)} CaloParticles")
    if "summary" in results:
        p = results["summary"]
        n0, n1 = p.n_fc0.mean(), p.n_fc1.mean()
        print(f"SimClusters per CaloParticle: Standard {n0:.2f}, FineCalo {n1:.2f} "
              f"(ratio {n1 / n0:.2f}x)")
        print(f"  same count in both configs: {100 * (p.n_fc0 == p.n_fc1).mean():.1f}% "
              f"of CaloParticles, FineCalo more: {100 * (p.n_fc1 > p.n_fc0).mean():.1f}%")
    if "surfaces" in results:
        labels = results["surfaces"]
        print("surfaces")
        for key in samples.KEYS:
            counts = labels[labels.config == key]["class"].value_counts()
            parts = ", ".join(f"{c} {counts.get(c, 0)} ({100 * counts.get(c, 0) / counts.sum():.1f}%)"
                              for c in classify_surfaces.CLASSES)
            print(f"  {samples.LABELS[key]:9s} {parts}")
    if results.get("decompositions") is not None:
        shares = results["decompositions"][0]["class"].value_counts(normalize=True)
        print("classes  " + ", ".join(f"{c} {100 * shares.get(c, 0):.1f}%"
                                      for c in shares.index))
    if results.get("compare"):
        print(f"compare  {len(results['compare'])} figures with both collections in "
              f"{Path(results['compare'][0]).parent}/")
    if "figures" in results or "surfaces" in results:
        n = len(results.get("figures", []))
        print(f"figures  {table_dir(sample.tree('fc0'), table)}/" + (f"  ({n} from the figures step)" if n else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", type=Path,
                    help="one nanoML file (its partner is found) or two, in either order")
    ap.add_argument("--table", choices=samples.TABLES,
                    help="default: boundary where both files carry it, else legacy")
    ap.add_argument("--tag", help="label for the output files; default: an existing "
                    "mapping's tag for these files, else the Standard file's stem")
    ap.add_argument("--steps", nargs="+", choices=STEPS, default=list(STEPS),
                    help="run only these steps (in the standard order); summary needs "
                    "the mappings from match, of this or an earlier run")
    ap.add_argument("--only", nargs="+", choices=make_figures.FIGURES, metavar="FIGURE",
                    default=list(make_figures.FIGURES),
                    help="figures step: only these, of " + ", ".join(make_figures.FIGURES))
    ap.add_argument("--excel", action="store_true",
                    help="also write the surface tables as a workbook next to the figures")
    ap.add_argument("--root", type=Path,
                    help="write to ROOT/outputs and ROOT/data/processed instead of the "
                    "repo's own, for a second production kept in its own folder")
    args = ap.parse_args()
    if len(args.files) > 2:
        ap.error("one or two files")

    if args.root:
        with output_paths.rooted(args.root):
            run(args, ap)
    else:
        run(args, ap)


def run(args, ap):
    sample = samples.load_sample(*args.files, tag=args.tag)
    table = args.table or sample.table
    if table not in sample.tables:
        ap.error(f"{table} table not in both files ({sample.tables})")
    suffix = "_bound" if table == "boundary" else ""
    out = table_dir(sample.tree("fc0"), table)
    out.mkdir(parents=True, exist_ok=True)
    report = out / f"report_{sample.tag}{suffix}.txt"

    steps = {
        "match": lambda: match_events.run(sample),
        "summary": lambda: paired_calopart_summary.run(sample, table),
        "figures": lambda: make_figures.run(sample, table, args.only),
        "surfaces": lambda: classify_surfaces.run(sample, table, excel=args.excel),
        "decompositions": lambda: classify_decompositions.run(sample, table),
        "compare": lambda: (make_figures.compare(sample, args.only) + classify_surfaces.compare(sample)
                            if {"boundary", "legacy"} <= set(sample.tables)
                            else print("Comparison figures need both SimCluster tables: skipped")),
    }
    results = {}
    with open(report, "w") as f, redirect_stdout(Tee(sys.stdout, f)):
        print(f"analyze.py {' '.join(str(p) for p in args.files)}")
        for i, step in enumerate(s for s in STEPS if s in args.steps):
            banner(f"Step {i + 1}: {step}")
            results[step] = steps[step]()
        headline(sample, table, results)
        print(f"report   {report}")


if __name__ == "__main__":
    main()
