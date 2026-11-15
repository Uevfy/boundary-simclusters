"""Compare each paired CaloParticle's decomposition across Standard and FineCalo, exactly.

The step after `match_events`: within a paired CaloParticle, the SimClusters of
the two configs are compared through the Geant4 track that opened each one
(`root_track`). Both configs simulate the same physics from the same seed and
differ only in what they store, so a track stored in both files carries the same
ID. This is measured here on every paired event (`id_agreement`) and the step
stops if it ever fails, so it has to hold for any new production too. Every
SimCluster counts, on whichever face it lies.

Where the two sets of root tracks differ, the parent links say how
(data/processed/ancestry_<tag>_<config>.parquet, written from RECO by
`write_ancestry`; nanoML has none). scripts/classify.py writes them and runs this. A FineCalo-only root is followed up its parent chain
to a Standard-only root. Each CaloParticle then has one class:

    same tracks         both configs record the same root tracks
    photon conversion   every replaced Standard root is a photon whose FineCalo
                        descendants come through its conversion
    other replacement   a Standard root is replaced by its FineCalo descendants,
                        at least one not through a photon's conversion: an
                        object lost in the gap
    no descendant       a Standard root has no FineCalo descendant, and nothing
                        is replaced
    mixed               a replacement and a Standard root with no descendant
    unlinked            a FineCalo-only root with no Standard-only ancestor
    different tracks    the root tracks differ and no ancestry is available

The legacy table keeps a crossing particle's crossing descendants too, so a
FineCalo-only root can also descend from a root both configs record, and Standard's
first-step rule can flag a pair below a photon that only FineCalo records. It
has its own classes, in `LEGACY_CLASSES`:

    same tracks         both configs record the same root tracks
    additions           every Standard root is in FineCalo too, which adds descendants
    replacement         a Standard root is replaced by its FineCalo descendants
    below finecalo      a Standard-only root descends from a FineCalo-only one
    unlinked            a FineCalo-only root with no Standard ancestor

Output: data/processed/decompositions_<tag>[_bound].parquet (one row per paired
CaloParticle) and decomposition_simclusters_<tag>[_bound].parquet (one row per
SimCluster: its role, the Standard root it descends from, where it was born).
"""

import numpy as np
import pandas as pd

from boundary_simclusters import output_paths, samples
from boundary_simclusters.boundary_geometry import label_crossings
from boundary_simclusters.reco_simtracks import check_alignment, read

KEYS = samples.KEYS
CLASSES = ("same tracks", "photon conversion", "other replacement", "no descendant", "mixed",
           "unlinked", "different tracks")
LEGACY_CLASSES = ("same tracks", "additions", "replacement", "below finecalo", "unlinked",
                  "no descendant", "mixed", "different tracks")
CONVERSION = 14                 # Geant4 process subtype of SimVertex.procType
PHOTON_ENERGY_TOLERANCE = 1e-5  # relative


def ancestry_path(tag, key):
    """Where config `key`'s parent links of the pair `tag` are written."""
    return output_paths.PROCESSED / f"ancestry_{tag.removesuffix('legacy')}_{key}.parquet"


def ancestry_table(tracks, vertices):
    """One row per Geant4 track: its parent, species, creating process and birth
    point. `event` is the row of the nanoML file the RECO file lines up with."""
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


def write_ancestry(sample, reco_paths):
    """Read each RECO file once, tie it to the config whose nanoML rows it lines
    up with (`check_alignment`; RECO file names say nothing reliable), and write
    that config's parent links. Returns {config: written path}."""
    output_paths.PROCESSED.mkdir(parents=True, exist_ok=True)
    nano = {k: sample.tree(k).arrays(["SimTrack_trackId", "SimTrack_pdgId"], library="np")
            for k in KEYS}
    written = {}
    for path in reco_paths:
        print(f"reading {path.name} (a few minutes)")
        tracks, vertices = read(path)
        for key in (k for k in KEYS if k not in written):
            try:
                check_alignment(tracks, nano[key])
            except SystemExit:
                continue
            out = ancestry_path(sample.tag, key)
            ancestry_table(tracks, vertices).to_parquet(out, index=False)
            print(f"{samples.LABELS[key]}: {path.name} -> {out.name}")
            written[key] = out
            break
    missing = [samples.LABELS[k] for k in KEYS if k not in written]
    if missing:
        raise SystemExit(f"no RECO file given lines up with the {' or '.join(missing)} nanoML file")
    return written


def read_ancestry(tag, key):
    """{event: {track: (parent, process, pdg, birth |z|)}} for one config, or None."""
    tag = tag.removesuffix("legacy")     # one ancestry per pair, whichever table is read
    for folder in (output_paths.PROCESSED, output_paths.REPO / "data" / "processed"):
        path = folder / f"ancestry_{tag}_{key}.parquet"
        if path.exists():
            d = pd.read_parquet(path)
            return {event: dict(zip(g.track.tolist(),
                                    zip(g.parent.tolist(), g.process.tolist(), g.pdg.tolist(),
                                        g.birth_abs_z.tolist())))
                    for event, g in d.groupby("event")}
    return None


def id_agreement(sample, caloparts):
    """(same, different): tracks stored under one ID in both files of a paired event,
    with identical or differing species and production pt, eta, phi."""
    names = [f"SimTrack_{b}" for b in ("trackId", "pdgId", "pt", "eta", "phi")]
    arrays = {k: sample.tree(k).arrays(names, library="np") for k in KEYS}

    def tracks(key, row):
        return dict(zip(arrays[key][names[0]][row].tolist(),
                        zip(*(arrays[key][n][row].tolist() for n in names[1:]))))

    same = different = 0
    for fc0_row, fc1_row in caloparts[["fc0_row", "fc1_row"]].drop_duplicates().itertuples(index=False):
        a, b = tracks("fc0", fc0_row), tracks("fc1", fc1_row)
        for track in a.keys() & b.keys():
            if a[track] == b[track]:
                same += 1
            else:
                different += 1
    return same, different


def climb(parents, track, targets):
    """(nearest ancestor of `track` in `targets`, the child just below it), or (None, None)."""
    child = track
    while True:
        link = parents.get(child)
        if link is None or link[0] == -1:
            return None, None
        if link[0] in targets:
            return link[0], child
        child = link[0]


def classify(std_only, fine_only, std, fine, parents):
    """(class, {root: role}, {FineCalo root: (ancestor, conversion?)}) of one
    CaloParticle with differing root tracks."""
    if parents is None:
        return "different tracks", {}, {}
    below = {f: climb(parents, f, std_only) for f in fine_only}
    links, used = {}, {}
    for f, (ancestor, child) in below.items():
        if ancestor is None:
            continue
        is_conversion = (std[ancestor] == 22 and parents[child][1] == CONVERSION)
        links[f] = (ancestor, is_conversion)
        used[ancestor] = used.get(ancestor, True) and is_conversion
    lonely = std_only - set(used)
    if len(links) < len(fine_only):
        cls = "unlinked"
    elif lonely and used:
        cls = "mixed"
    elif lonely:
        cls = "no descendant"
    elif all(used.values()):
        cls = "photon conversion"
    else:
        cls = "other replacement"
    roles = {s: "replaced" if s in used else "no descendant" for s in std_only}
    roles.update({f: "descendant" if f in links else "unlinked" for f in fine_only})
    return cls, roles, links


def climbs_to(parents, track, targets):
    """Does `track`'s parent chain reach `targets`?"""
    return climb(parents, track, targets)[0] is not None


def classify_legacy(std_only, fine_only, std, fine, parents_std, parents_fine):
    """The same, on the legacy table, where a root can descend from a shared root too."""
    if parents_fine is None:
        return "different tracks", {}, {}
    standard = set(std)
    links = {f: climb(parents_fine, f, standard)[0] for f in fine_only}
    used = {a for a in links.values() if a is not None}
    below = {s for s in std_only if climbs_to(parents_std, s, fine_only) or climbs_to(parents_fine, s, fine_only)}
    lonely = std_only - used - below
    roles = {s: "below finecalo" if s in below else "replaced" if s in used else "no descendant" for s in std_only}
    roles.update({f: "unlinked" if a is None else "descendant" if a in std_only else "addition"
                  for f, a in links.items()})
    if below:
        cls = "below finecalo"
    elif None in links.values():
        cls = "unlinked"
    elif lonely:
        cls = "mixed" if fine_only else "no descendant"
    elif std_only:
        cls = "replacement"
    else:
        cls = "additions"
    return cls, roles, {f: (a, False) for f, a in links.items() if a is not None}


def compare(sample, tag, table):
    """(one row per CaloParticle, one row per SimCluster)."""
    _, caloparts = sample.mappings()
    frames = {k: label_crossings(samples.simcluster_frame(sample.tree(k), table, k), k) for k in KEYS}
    groups = {k: dict(tuple(frames[k].groupby(["event", "cp"]))) for k in KEYS}
    empty = {k: frames[k].iloc[:0] for k in KEYS}
    ancestry = read_ancestry(tag, "fc1")
    ancestry_std = read_ancestry(tag, "fc0") if table == "legacy" else None
    if ancestry is None:
        print(f"WARNING: no ancestry parquet for {tag}, so every differing CaloParticle is "
              "'different tracks'. Run scripts/classify.py with the RECO files.")
    rows, clusters = [], []
    for n, cp in enumerate(caloparts.itertuples()):
        a = groups["fc0"].get((cp.fc0_row, cp.fc0_cp), empty["fc0"])
        b = groups["fc1"].get((cp.fc1_row, cp.fc1_cp), empty["fc1"])
        std = dict(zip(a.root_track.tolist(), a.pdg.tolist()))
        fine = dict(zip(b.root_track.tolist(), b.pdg.tolist()))
        for track in std.keys() & fine.keys():
            if std[track] != fine[track]:
                raise SystemExit(f"{sample.tag}: track {track} changes species across configs")
        std_only, fine_only = set(std) - set(fine), set(fine) - set(std)
        if std_only or fine_only:
            parents = None if ancestry is None else ancestry.get(cp.fc1_row, {})
            if table == "legacy":
                above = {} if ancestry_std is None else ancestry_std.get(cp.fc0_row, {})
                cls, roles, links = classify_legacy(std_only, fine_only, std, fine, above, parents)
            else:
                cls, roles, links = classify(std_only, fine_only, std, fine, parents)
        else:
            cls, roles, links = "same tracks", {}, {}
        rows.append(dict(calopart=n, **{"class": cls}, n_std=len(a), n_fine=len(b),
                         n_shared=len(set(std) & set(fine)), n_std_only=len(std_only),
                         n_fine_only=len(fine_only),
                         n_side_entries=int(((b["group"] == "inner side") & b["entering"]).sum())))
        for key, d in (("fc0", a), ("fc1", b)):
            clusters.append(d.assign(
                config=key, calopart=n,
                role=[roles.get(t, "shared") for t in d.root_track.tolist()],
                ancestor=[links.get(t, (-1, False))[0] for t in d.root_track.tolist()],
                through_conversion=[links.get(t, (-1, False))[1] for t in d.root_track.tolist()]))
    own = samples.caloparticle_frame(sample.tree("fc0")).set_index(["event", "cp"]).loc[
        list(zip(caloparts.fc0_row, caloparts.fc0_cp))]
    result = pd.DataFrame(rows).assign(
        fc0_row=caloparts.fc0_row.to_numpy(), fc0_cp=caloparts.fc0_cp.to_numpy(),
        fc1_row=caloparts.fc1_row.to_numpy(), fc1_cp=caloparts.fc1_cp.to_numpy(),
        pdg=own.pdg.to_numpy().astype(int), abs_eta=np.abs(own.eta.to_numpy()),
        energy=own.energy.to_numpy())
    columns = ["config", "calopart", "event", "cp", "root_track", "pdg", "ekin", "role", "ancestor",
               "through_conversion", "root_is_primary", "group", "entering", "abs_z", "r"]
    return result, pd.concat(clusters, ignore_index=True)[columns]


def summarise(result, clusters, ids):
    """The numbers `pipeline.numbers` carries."""
    a, b = (clusters[clusters.config == k] for k in KEYS)
    shared = a[a.role == "shared"].merge(b[b.role == "shared"], on=["calopart", "root_track"],
                                         suffixes=("_std", "_fine"))
    photons = shared[shared.pdg_std == 22]
    charged = shared[(shared.pdg_std != 22) & (shared.ekin_std > 0)]
    differs = np.abs(photons.ekin_fine - photons.ekin_std) > PHOTON_ENERGY_TOLERANCE * photons.ekin_std
    return {
        "class": {str(k): int(v) for k, v in sorted(result["class"].value_counts().items())},
        "tracks_same_id": ids[0],
        "tracks_id_disagree": ids[1],
        "shared_simclusters": len(shared),
        "shared_photons": len(photons),
        "shared_photon_energy_differs": int(differs.sum()),
        "shared_charged_ekin_change_median": float(
            np.median((charged.ekin_fine - charged.ekin_std) / charged.ekin_std)) if len(charged) else 0.0,
        "standard_only": int((a.role != "shared").sum()),
        "standard_replaced": int((a.role == "replaced").sum()),
        "standard_below_finecalo": int((a.role == "below finecalo").sum()),
        "standard_no_descendant": int((a.role == "no descendant").sum()),
        "finecalo_only": int((b.role != "shared").sum()),
        "finecalo_descendants": int((b.role == "descendant").sum()),
        "finecalo_additions": int((b.role == "addition").sum()),
        "finecalo_unlinked": int((b.role == "unlinked").sum()),
    }


def report(result, numbers, classes=CLASSES):
    print(f"\nTracks stored under one ID in both files: {numbers['tracks_same_id']} identical, "
          f"{numbers['tracks_id_disagree']} different")
    print("\nCaloParticles by class")
    print("----------------------")
    for cls in classes:
        n = int((result["class"] == cls).sum())
        print(f"  {cls:18s} {n:5d}  {100 * n / len(result):5.2f}%")
    print(f"\nFineCalo-only SimClusters {numbers['finecalo_only']}: "
          f"{numbers['finecalo_descendants']} descend from a Standard-only one, "
          f"{numbers['finecalo_unlinked']} do not")
    print(f"Standard-only SimClusters {numbers['standard_only']}: "
          f"{numbers['standard_replaced']} with a FineCalo descendant, "
          f"{numbers['standard_no_descendant']} without")
    print(f"Shared photons {numbers['shared_photons']}: energy differs in "
          f"{numbers['shared_photon_energy_differs']}; shared charged, kinetic energy from "
          f"Standard's surface to FineCalo's, median "
          f"{100 * numbers['shared_charged_ekin_change_median']:+.1f}%  (the gap, not the mechanism)")


def run(sample, table=None, tag=None):
    """Classify `sample`'s CaloParticles, print the summary, write both tables,
    and return (per CaloParticle, numbers)."""
    table_name = table or sample.table
    tag = tag or sample.tag
    _, caloparts = sample.mappings()
    ids = id_agreement(sample, caloparts)
    if ids[1]:
        raise SystemExit(f"{sample.tag}: {ids[1]} tracks carry one ID but differ across the configs; "
                         "track IDs cannot be compared on this pair")
    result, clusters = compare(sample, tag, table_name)
    numbers = summarise(result, clusters, ids)
    report(result, numbers, LEGACY_CLASSES if table_name == "legacy" else CLASSES)
    suffix = "_bound" if table_name == "boundary" else ""
    result.to_parquet(output_paths.PROCESSED / f"decompositions_{tag}{suffix}.parquet", index=False)
    clusters.to_parquet(output_paths.PROCESSED / f"decomposition_simclusters_{tag}{suffix}.parquet", index=False)
    return result, numbers
