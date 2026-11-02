"""Test the boundary-SimCluster rule against Geant4's full ancestry.

`CaloTruthAccumulator` (stock, `dgaytanv/cmssw` `pepr_CMSSW_20_0_0_pre1`) opens
a boundary SimCluster at the first `crossedBoundary()` track on each branch of a
CaloParticle's decay tree (and only if the subtree has simhits), then absorbs
everything below it without testing the selector again. The prediction:

    a crossing track gets its own boundary SimCluster
        <=> no ancestor of it crossed (and its subtree left hits)

and absorbed crossings keep their hits in the ancestor's cluster but not their
boundary energy. That is the boundary-energy deficit (SimTrack minus
SimCluster): crossings absorbed that were born *before*
the boundary (radiated photons, mostly) carry energy recorded nowhere, while
those born *inside* the volume (conversion pairs of a photon that had already
crossed) carry energy their ancestor's cluster already counts.

Nano has no parent branches, so the ancestry comes from the RECO file's
`g4SimHits` collections via `reco_simtracks.py`. The nano file must be the
nanoML of that same RECO file: the rows are checked to line up (every nano
SimTrack found in the same RECO event with the same pdgId) before anything is
counted. The two files can be given in either order: which is RECO and which
nanoML is read from their branches, and whose boundary volume decides "born
inside" from the nanoML file's config (`samples.py`). It is the one script in
the repo that reads a RECO file.

Usage:
    python scripts/studies/check_boundary_nesting.py photon_fc0_20GeV_testRECO.root \
        bound_photon_fc0_20GeV_testNanoML.root
"""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import uproot

from boundary_simclusters import samples
from boundary_simclusters.boundary_geometry import radius, source_inside
from boundary_simclusters.reco_simtracks import (CROSSED_BOUNDARY, TRACKS, ancestry, check_alignment,
                                                 read)

LABELS = samples.LABELS
NANO_BRANCHES = ["SimTrack_trackId", "SimTrack_pdgId", "BoundarySimCluster_trackIdAtBoundary",
                 "SimCluster_trackIdAtBoundary"]


def sort_files(first, second):
    """(RECO, nanoML) from two paths in either order, told apart by their branches."""
    kinds = {}
    for path in (first, second):
        keys = set(uproot.open(path)["Events"].keys(recursive=False))
        if TRACKS in keys:
            kinds.setdefault("reco", path)
        elif "SimTrack_trackId" in keys:
            kinds.setdefault("nano", path)
    if set(kinds) != {"reco", "nano"}:
        raise SystemExit(f"need one RECO file (with {TRACKS.rstrip('.')}) and one nanoML file; "
                         f"got {first.name}, {second.name}")
    return kinds["reco"], kinds["nano"]


def kinetic_at_boundary(t):
    p4 = t["momAtBoundary"].astype(np.float64)
    e, p2 = p4[:, 3], (p4[:, :3] ** 2).sum(axis=1)
    return e - np.sqrt(np.clip(e ** 2 - p2, 0.0, None))


def classify(tracks, vertices, nano, config):
    """One row per crossing track: (event, trackId, pdgId, E_kin at boundary,
    has a crossing ancestor, ancestor's pdgId, born inside the volume, is a
    boundary-SimCluster root, has a legacy SimCluster)."""
    parents = ancestry(tracks, vertices)
    rows = []
    for i, (t, v, parent) in enumerate(zip(tracks, vertices, parents)):
        tid = t["trackId"].astype(int).tolist()
        crossed = dict(zip(tid, ((t["trackInfo"] & CROSSED_BOUNDARY) > 0).tolist()))
        pdg = dict(zip(tid, t["pdgId"].tolist()))
        vertex = v["pos"][t["ivert"]]
        inside = dict(zip(tid, source_inside(np.abs(vertex[:, 2]),
                                             radius(vertex[:, 0], vertex[:, 1]), config).tolist()))
        ekin = dict(zip(tid, kinetic_at_boundary(t).tolist()))
        roots = set(nano["BoundarySimCluster_trackIdAtBoundary"][i].astype(int).tolist())
        legacy = set(nano["SimCluster_trackIdAtBoundary"][i].astype(int).tolist())
        for x in tid:
            if not crossed[x]:
                if x in roots:
                    raise RuntimeError(f"event {i}: SimCluster root {x} did not cross")
                continue
            a = parent[x]
            while a != -1 and a in parent and not crossed[a]:
                a = parent[a]
            nested = a != -1 and a in parent
            rows.append((i, x, pdg[x], ekin[x], nested, pdg[a] if nested else 0,
                         inside[x], x in roots, x in legacy))
    return rows


def report(rows, n_events, config):
    print(f"\n{LABELS[config]}: {len(rows)} crossing tracks in {n_events} events\n")
    table = Counter((r[4], r[7]) for r in rows)
    print("  crossing track             gets a boundary SimCluster   does not")
    for nested, name in ((True, "has a crossing ancestor"), (False, "no ancestor crossed")):
        print(f"  {name:26s} {table[(nested, True)]:>27d}   {table[(nested, False)]:>8d}")
    verdict = "holds" if table[(True, True)] == 0 else "VIOLATED"
    print(f"\n  rule (no SimCluster below a crossing): {verdict}")

    free = [r for r in rows if not r[4] and not r[7]]
    if free:
        print(f"  first crossings without a SimCluster: {len(free)}, max E_kin "
              f"{max(r[3] for r in free):.3f} GeV, with hits of their own (legacy SimCluster): "
              f"{sum(r[8] for r in free)}  -- the selector's simhits condition")

    def per_event(sel):
        """Summed kinetic energy of the selected tracks in each event (0 if none)."""
        sums = np.zeros(n_events)
        for r in rows:
            if sel(r):
                sums[r[0]] += r[3]
        return sums

    recorded = per_event(lambda r: r[7])
    before = per_event(lambda r: r[4] and not r[6])
    inside = per_event(lambda r: r[4] and r[6])
    # Means add up line by line; medians do not, so the total's median is its own.
    print("\n  kinetic energy at the boundary, per event      mean      median")
    for name, e in (("recorded in the boundary table (roots)", recorded),
                    ("absorbed, born before the boundary", before),
                    ("= total crossing the boundary", recorded + before),
                    ("absorbed, born inside the volume (not added)", inside)):
        print(f"    {name:44s} {e.mean():7.3f}   {np.median(e):7.3f} GeV")

    upstream = [r for r in rows if r[4] and not r[6]]
    species = Counter(r[2] for r in upstream).most_common(4)
    ancestors = Counter(r[5] for r in upstream).most_common(4)
    print(f"\n  absorbed before the boundary: {len(upstream)} tracks, pdgId {species}")
    print(f"    their first crossing ancestor: pdgId {ancestors}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs=2, type=Path,
                    help="a RECO file and the nanoML file made from it, in either order")
    args = ap.parse_args()

    reco, nano_path = sort_files(*args.files)
    tree = samples.events(nano_path)
    config = samples.detect_config(tree)
    tracks, vertices = read(reco)
    nano = tree.arrays(NANO_BRANCHES, library="np")
    check_alignment(tracks, nano)
    print(f"{reco.name}: {sum(len(t['trackId']) for t in tracks)} SimTracks, "
          f"{sum(len(v['vtxId']) for v in vertices)} SimVertices; rows line up with {nano_path.name}")
    count = lambda branch: sum(len(x) for x in nano[branch])
    print(f"SimClusters in {nano_path.name}: boundary table "
          f"{count('BoundarySimCluster_trackIdAtBoundary')}, legacy table "
          f"{count('SimCluster_trackIdAtBoundary')}")
    report(classify(tracks, vertices, nano, config), len(tracks), config)


if __name__ == "__main__":
    main()
