"""Classify each CaloParticle's decomposition: nanoML and RECO files in, classes out.

    python scripts/classify.py photon_20GeV_fc0_NanoML.root photon_20GeV_fc1_NanoML.root \\
        photon_20GeV_fc0_RECO.root photon_20GeV_fc1_RECO.root
    python scripts/classify.py photon_20GeV_fc0_NanoML.root      # ancestry already written
    python scripts/classify.py a.root b.root c.root d.root --table legacy

The second half of the standard run; scripts/analyze.py is the first and needs
nanoML only. Within each paired CaloParticle, the SimClusters of Standard and
FineCalo are compared through the Geant4 track that opened each, and where they
differ, the parent links from RECO say how: one class per CaloParticle
(`pipeline.classify_decompositions`).

Give the files in any order: RECO and nanoML are told apart by their contents,
and each RECO file is tied to the nanoML file it was made from by its rows.
Give one nanoML file and its partner is found in the same folder.

    ancestry   nanoML has no parent links, so they are read from RECO once
               (minutes per file) into data/processed/ancestry_<tag>_<fc0|fc1>.parquet
               and reused afterwards: later runs need no RECO. Delete those
               files to read RECO again.
    match      the event and CaloParticle pairing of analyze.py's match step,
               run here if it has not been yet
    classes    data/processed/decompositions_<tag>[_bound].parquet (one row per
               CaloParticle) and decomposition_simclusters_<tag>[_bound].parquet

Everything printed is also written to report_<tag>[_bound]_classes.txt, beside
analyze.py's report. `--tag`, `--table` and `--root` as in scripts/analyze.py.
"""

import argparse
import sys
from contextlib import redirect_stdout
from pathlib import Path

import uproot

from boundary_simclusters import output_paths, samples
from boundary_simclusters.output_paths import table_dir
from boundary_simclusters.pipeline import classify_decompositions, match_events
from boundary_simclusters.reco_simtracks import TRACKS


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


def is_reco(path):
    """A RECO file carries Geant4's SimTracks; a nanoML file does not."""
    return TRACKS in uproot.open(path)["Events"].keys(recursive=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", type=Path,
                    help="one or two nanoML files and, the first time, both RECO files; any order")
    ap.add_argument("--table", choices=samples.TABLES,
                    help="default: boundary where both files carry it, else legacy")
    ap.add_argument("--tag", help="as for analyze.py; default <particle>_<energy>, e.g. photon_20GeV")
    ap.add_argument("--root", type=Path,
                    help="write to ROOT/outputs and ROOT/data/processed, as for analyze.py")
    args = ap.parse_args()
    if args.root:
        with output_paths.rooted(args.root):
            run(args, ap)
    else:
        run(args, ap)


def run(args, ap):
    reco = [f for f in args.files if is_reco(f)]
    nano = [f for f in args.files if f not in reco]
    if not 1 <= len(nano) <= 2 or len(reco) > 2:
        ap.error(f"need one or two nanoML files and at most two RECO files "
                 f"(got {len(nano)} nanoML, {len(reco)} RECO)")
    sample = samples.load_sample(*nano, tag=args.tag)
    table = args.table or sample.table
    if table not in sample.tables:
        ap.error(f"{table} table not in both files ({sample.tables})")
    suffix = "_bound" if table == "boundary" else ""
    out = table_dir(sample.tree("fc0"), table)
    out.mkdir(parents=True, exist_ok=True)
    report = out / f"report_{sample.tag}{suffix}_classes.txt"

    with open(report, "w") as f, redirect_stdout(Tee(sys.stdout, f)):
        print(f"classify.py {' '.join(str(p) for p in args.files)}")
        print(f"sample {sample.tag}: {sample.particle} {sample.energy}, {table} table")

        banner("Step 1: ancestry")
        paths = {k: classify_decompositions.ancestry_path(sample.tag, k) for k in samples.KEYS}
        if all(p.exists() for p in paths.values()):
            print("using " + ", ".join(p.name for p in paths.values())
                  + (" (RECO files not read; delete these to read them again)" if reco else ""))
        elif len(reco) == 2:
            classify_decompositions.write_ancestry(sample, reco)
        else:
            raise SystemExit(f"no ancestry for {sample.tag} yet: give both RECO files as well")

        banner("Step 2: match")
        mapping = output_paths.PROCESSED / f"event_matching_{sample.tag}.parquet"
        if mapping.exists():
            sample.mappings()            # stops if the mapping belongs to other files
            print(f"using the mappings of an earlier match run ({mapping.name})")
        else:
            match_events.run(sample)

        banner("Step 3: classes")
        result, _ = classify_decompositions.run(sample, table)

        banner("Summary")
        shares = result["class"].value_counts(normalize=True)
        print(f"sample   {sample.tag}: {sample.particle} {sample.energy}, {table} table, "
              f"{len(result)} CaloParticles")
        print("classes  " + ", ".join(f"{c} {100 * shares.get(c, 0):.1f}%" for c in shares.index))
        print(f"report   {report}")


if __name__ == "__main__":
    main()
