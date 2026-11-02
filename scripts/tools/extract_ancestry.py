"""Write each config's parent links from RECO to a parquet, for the decomposition classes.

    python scripts/tools/extract_ancestry.py bound_photon_fc0_20GeV_testNanoML.root \
        bound_photon_fc1_20GeV_testNanoML.root photon_fc0_20GeV_testRECO.root \
        photon_fc1_20GeV_testRECO.root

Takes a pair's two nanoML files and two RECO files, in any order; each RECO is
tied to the nanoML whose rows it lines up with (`reco_simtracks.check_alignment`),
since RECO file names say nothing reliable. NanoML has no parent links, and
reading RECO takes minutes per file, so the classification reads this instead:
data/processed/ancestry_<tag>_<fc0|fc1>.parquet, one row per stored track, with
`event` the row of the nanoML file. `--root` as in scripts/analyze.py.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

from boundary_simclusters import output_paths, samples
from boundary_simclusters.reco_simtracks import TRACKS, check_alignment, read


def table(tracks, vertices):
    frames = []
    for event, (t, v) in enumerate(zip(tracks, vertices)):
        vertex = t["ivert"]
        frames.append(pd.DataFrame({
            "event": event,
            "track": t["trackId"].astype(np.int64),
            "parent": v["parentTrackId"][vertex].astype(np.int64),
            "pdg": t["pdgId"].astype(np.int32),
            "process": v["procType"][vertex].astype(np.int32),
            "birth_abs_z": np.abs(v["pos"][vertex, 2]).astype(float),
            "birth_r": np.hypot(v["pos"][vertex, 0], v["pos"][vertex, 1]).astype(float)}))
    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs=4, type=Path, help="both nanoML and both RECO files, any order")
    ap.add_argument("--tag")
    ap.add_argument("--root", type=Path, help="write under this tree instead of the repo's")
    args = ap.parse_args()
    reco = [f for f in args.files if TRACKS in uproot.open(f)["Events"].keys(recursive=False)]
    nano = [f for f in args.files if f not in reco]
    if len(reco) != 2 or len(nano) != 2:
        raise SystemExit("need two RECO and two nanoML files")
    sample = samples.load_sample(*nano, tag=args.tag)
    processed = (args.root / "data" / "processed") if args.root else output_paths.PROCESSED
    for key in samples.KEYS:
        arrays = sample.tree(key).arrays(["SimTrack_trackId", "SimTrack_pdgId"], library="np")
        for path in reco:
            tracks, vertices = read(path)
            try:
                check_alignment(tracks, arrays)
            except SystemExit:
                continue
            out = processed / f"ancestry_{sample.tag}_{key}.parquet"
            table(tracks, vertices).to_parquet(out, index=False)
            print(f"{sample.tag} {samples.LABELS[key]}: {path.name} -> {out.name}")
            reco.remove(path)
            break
        else:
            raise SystemExit(f"no RECO file lines up with {sample.files[key].name}")


if __name__ == "__main__":
    main()
