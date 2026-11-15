"""What a nanoML file is, which file is its partner, and its SimClusters as one table.

Every fact about a sample is read from the file's contents, never from its name
(the numeric prefixes contradict the pdgId, README):

    config     Standard (fc0) or FineCalo (fc1): whichever volume's main face the
               SimClusters land on -- CALO's tracker-cavity step at |z| = 293.5 cm
               or HGCal's front face at 321.05 cm
    particle   GenPart_pdgId, energy CaloPart_energy (`output_paths`)
    tables     which SimCluster collections the file carries; the boundary table
               (`BoundarySimCluster_*`) only exists in the `bound_*` production
    partner    the file of the other config that holds the same events: the one
               whose GEN keys overlap

so a new file needs no edit anywhere before it can be analysed.

`simcluster_frame` is the one per-SimCluster reader: kinematics, the envelope face
it landed on, kinetic energy (`boundaryEnergy` minus rest mass) and the ancestry the
nano file can support (`ancestry`). Mappings between the two files of a pair come
from the match step (`pipeline.match_events`) and are found by the file names they record.

Imported by the analysis scripts. Run on its own, it says what files are:

    python -m boundary_simclusters.samples bound_photon_fc0_20GeV_testNanoML.root
    python -m boundary_simclusters.samples a.root b.root
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

from boundary_simclusters.boundary_geometry import MAIN_FACE, classify_faces, kinetic_energy, radius  # noqa: E402
from boundary_simclusters.matching import genpart_signature, read_mapping  # noqa: E402
from boundary_simclusters import output_paths
from boundary_simclusters.output_paths import energy_label, particle_name  # noqa: E402

KEYS = ("fc0", "fc1")
LABELS = {"fc0": "Standard", "fc1": "FineCalo"}
TABLES = {"boundary": "BoundarySimCluster", "legacy": "SimCluster"}

# Events read to decide a file's config. Every SimCluster of a local file lands
# on its own config's main face far more often than on the other's (>= 94% vs 0%),
# so a few hundred events decide it with no ambiguity.
CONFIG_EVENTS = 200
MIN_MAIN_FACE_SHARE = 0.5

# Share of a file's GEN keys another file must contain to be its partner. True
# partners share 97-100%. Samples of different seeds can share most of
# their events too (~88% between two 200 GeV e⁻ seeds), so the partner is the best overlap,
# and overlaps closer than OVERLAP_TIE count as a tie.
MIN_GEN_OVERLAP = 0.5
OVERLAP_TIE = 0.005


# --- what a file is ---------------------------------------------------------


def events(path):
    return uproot.open(Path(path))["Events"]


def has_table(tree, table):
    return f"{TABLES[table]}_pdgId" in tree.keys()


def detect_config(tree, n_events=CONFIG_EVENTS):
    """"fc0" or "fc1": the config whose main face the SimClusters land on.

    Read from the legacy table, which every nanoML file carries."""
    arr = tree.arrays([f"SimCluster_impactPoint_{a}" for a in "xyz"],
                      entry_stop=n_events, library="np")
    x, y, z = (np.concatenate(arr[f"SimCluster_impactPoint_{a}"]).astype(float)
               for a in "xyz")
    share = {}
    for key in KEYS:
        face, _ = classify_faces(np.abs(z), radius(x, y), key)
        share[key] = float(np.mean(face == MAIN_FACE[key])) if len(face) else 0.0
    best = max(share, key=share.get)
    if share[best] < MIN_MAIN_FACE_SHARE:
        raise ValueError(f"no config fits: main-face share {share}")
    return best


def describe(path):
    """Config, particle, energy and SimCluster tables of one nanoML file."""
    tree = events(path)
    return {
        "path": Path(path),
        "config": detect_config(tree),
        "particle": particle_name(tree),
        "energy": energy_label(tree),
        "tables": [t for t in TABLES if has_table(tree, t)],
        "events": tree.num_entries,
    }


def gen_overlap(tree_a, tree_b):
    """Share of `tree_a`'s GEN keys that also occur in `tree_b`."""
    a, b = genpart_signature(tree_a), set(genpart_signature(tree_b))
    return sum(key in b for key in a) / len(a) if a else 0.0


def find_partner(path, search_dir=None, candidates=None):
    """The file of the other config that holds the same events as `path`.

    Looks at `candidates` if given, else at every nanoML file in `search_dir`
    (default: the file's own folder). Of the files of the other config that share its events, keeps those
    with the highest GEN overlap, then those exporting the same SimCluster
    tables (both files of one production do). Raises if none is left or several
    are -- one gun seed can be produced more than once, and choosing between
    those is not something to guess.
    """
    path = Path(path).resolve()
    me = describe(path)
    tree = events(path)
    matches = []
    pool = candidates if candidates is not None else Path(search_dir or path.parent).glob("*.root")
    for other in sorted(Path(c) for c in pool):
        if other.resolve() == path:
            continue
        try:
            them = describe(other)
        except (KeyError, ValueError, OSError):
            continue                      # not a nanoML file (RECO, GSD, ...)
        if them["config"] == me["config"] or (them["particle"], them["energy"]) != (
                me["particle"], me["energy"]):
            continue
        overlap = gen_overlap(tree, events(other))
        if overlap >= MIN_GEN_OVERLAP:
            matches.append((overlap, them["tables"] == me["tables"], other))
    if not matches:
        where_ = "the given files" if candidates is not None else path.parent
        raise SystemExit(f"no partner for {path.name} in {where_}")
    best = max(o for o, _, _ in matches)
    left = [m for m in matches if m[0] >= best - OVERLAP_TIE]
    if any(same for _, same, _ in left):
        left = [m for m in left if m[1]]
    if len(left) > 1:
        names = ", ".join(f"{m.name} ({o:.1%})" for o, _, m in left)
        raise SystemExit(f"{path.name} has several possible partners: {names}; "
                         f"name the partner explicitly")
    return left[0][2]


# --- a matched pair ---------------------------------------------------------


@dataclass
class Sample:
    """One fc0/fc1 pair, every attribute read from the files."""
    tag: str
    files: dict          # {"fc0": Path, "fc1": Path}
    particle: str
    energy: str
    tables: list         # tables both files carry, boundary first

    @property
    def table(self):
        """The table the standard run reads: boundary where it exists."""
        return self.tables[0]

    def tree(self, key):
        return events(self.files[key])

    def mappings(self):
        """(event mapping, CaloParticle mapping) from the match step."""
        args = (self.tag, self.files["fc0"].name, self.files["fc1"].name)
        try:
            return (read_mapping("event_matching", *args),
                    read_mapping("calopart_matching", *args))
        except FileNotFoundError:
            raise SystemExit(
                f"no mapping for {self.tag}; run: python scripts/analyze.py "
                f"{self.files['fc0'].name} {self.files['fc1'].name} "
                f"--tag {self.tag} --steps match") from None


def existing_tag(fc0, fc1):
    """The tag the match step already used for this pair of files, if any."""
    for mapping in sorted(output_paths.PROCESSED.glob("event_matching_*.parquet")):
        files = pd.read_parquet(mapping, columns=["fc0_file", "fc1_file"]).iloc[0]
        if (files.fc0_file, files.fc1_file) == (Path(fc0).name, Path(fc1).name):
            return mapping.stem.removeprefix("event_matching_")
    return None


def load_sample(first, second=None, tag=None):
    """The pair `first` belongs to, in either order; its partner is found if not given.

    The tag is, in order: the one given, the one an existing mapping of these
    two files was written under, or the fc0 file's stem. It is a label only.
    """
    second = second or find_partner(first)
    found = {d["config"]: d for d in (describe(first), describe(second))}
    if set(found) != set(KEYS):
        raise SystemExit(f"{Path(first).name} and {Path(second).name} are both "
                         f"{LABELS[next(iter(found))]}")
    fc0, fc1 = found["fc0"], found["fc1"]
    for field in ("particle", "energy"):
        if fc0[field] != fc1[field]:
            raise SystemExit(f"not one sample: {field} {fc0[field]} vs {fc1[field]}")
    if gen_overlap(events(fc0["path"]), events(fc1["path"])) < MIN_GEN_OVERLAP:
        raise SystemExit(f"{fc0['path'].name} and {fc1['path'].name} do not hold "
                         f"the same events (GEN keys do not overlap)")
    tag = tag or existing_tag(fc0["path"], fc1["path"]) or fc0["path"].stem
    tables = [t for t in TABLES if t in fc0["tables"] and t in fc1["tables"]]
    return Sample(tag=tag, files={"fc0": fc0["path"], "fc1": fc1["path"]},
                  particle=fc0["particle"], energy=fc0["energy"], tables=tables)


# --- SimClusters ------------------------------------------------------------


def where(face, key):
    """Coarse placement: main face, another face of the envelope, or neither."""
    return np.where(face == MAIN_FACE[key], "main face",
                    np.where(np.isin(face, ("interior", "exterior")), face,
                             "other face"))


def relation(root_pdg, calopart_pdg, is_primary):
    """How the root track's species relates to its CaloParticle's.

    A *species pattern*, not a measured chain: "e± from γ" says a photon's
    CaloParticle is represented by an electron or positron at the boundary,
    which is what a conversion looks like, without the parent link that would
    prove it (docs/data_handling.md).
    """
    root, cp = np.abs(root_pdg), np.abs(calopart_pdg)
    return np.select(
        [is_primary,
         (cp == 22) & (root == 11),
         (cp == 11) & (root == 22),
         root == cp],
        ["primary", "e± from γ", "γ from e±", "same species, secondary"],
        default="other")


def ancestry(root_track, root_pdg, calopart_track, calopart_pdg, face):
    """What the nano file can say about how a SimCluster descends from its CaloParticle.

    root_track       Geant4 track that crossed and opened the SimCluster
    root_is_primary  that track is the CaloParticle's own primary
    relation         species pattern of root vs CaloParticle (`relation`)
    born_inside      flagged at its vertex inside the volume, on no face

    The chain between primary and root is not in nanoML: it has no parent links
    (docs/data_handling.md; production and decay vertices exist only for the
    tracks with a TrackingParticle, `TrackingPart_Vtx_*` / `_DecayVtx_*`), and
    `SimTrack_*SimClusterIdx` points
    only at each SimCluster's root track, never at what was folded into it.
    Parent branches added to nanoML_cfg.py would extend this function with the
    chain (depth, intermediate species); every column above stays valid.
    """
    is_primary = root_track == calopart_track
    return {
        "root_track": root_track,
        "root_is_primary": is_primary,
        "relation": relation(root_pdg, calopart_pdg, is_primary),
        "born_inside": face == "interior",
    }


def simcluster_frame(tree, table, key):
    """One row per SimCluster of `table` in a file of config `key`.

    With the kinematics the figure scripts read and the ancestry columns.
    """
    p = TABLES[table]
    jag = lambda b: tree[b].array(library="np")
    counts = np.array([len(x) for x in jag(f"{p}_pdgId")])
    flat = lambda b: np.concatenate(jag(f"{p}_{b}"))
    event = np.repeat(np.arange(len(counts)), counts)
    cp = flat("CaloPartIdx").astype(int)

    x, y, z = (flat(f"impactPoint_{a}").astype(float) for a in "xyz")
    r = radius(x, y)
    face, _ = classify_faces(np.abs(z), r, key)
    energy = flat("boundaryEnergy").astype(float)
    pmag = flat("boundaryPmag").astype(float)
    pdg = flat("pdgId")

    cp_track, cp_pdg = (jag(f"CaloPart_{b}") for b in ("trackId", "pdgId"))
    own = lambda per_event: np.array([per_event[e][c] for e, c in zip(event, cp)])

    frame = pd.DataFrame({
        "event": event,
        "cp": cp,
        "pdg": pdg,
        "x": x, "y": y, "z": z, "abs_z": np.abs(z), "r": r,
        "eta": flat("impactPoint_eta").astype(float),
        "phi": flat("impactPoint_phi").astype(float),
        "energy": energy,
        "pmag": pmag,
        "pt": flat("boundaryPt").astype(float),
        "ekin": kinetic_energy(energy, pmag),
        "nhits": flat("nHits"),
        "rec": flat("recEnergy").astype(float),
        "face": face,
        "where": where(face, key),
    })
    links = ancestry(flat("trackIdAtBoundary"), pdg, own(cp_track), own(cp_pdg), face)
    return frame.assign(**links)


def simtracks(tree):
    """Every SimTrack that crossed, with its crossing point and momentum there.
    No energy branch exists at track level, so no kinetic energy."""
    jag = tree["SimTrack_crossedBoundary"].array(library="np")
    flat = lambda b: np.concatenate(tree[b].array(library="np"))
    keep = np.concatenate(jag).astype(bool)
    event = np.repeat(np.arange(len(jag)), [len(e) for e in jag])[keep]
    x, y, z = (flat(f"SimTrack_boundaryPos_{a}")[keep].astype(float) for a in "xyz")
    r = radius(x, y)
    return pd.DataFrame({
        "event": event,
        "x": x, "y": y, "z": z, "abs_z": np.abs(z), "r": r,
        "pt": flat("SimTrack_boundaryMomentum_pt")[keep].astype(float),
        "eta": flat("SimTrack_boundaryMomentum_eta")[keep].astype(float),
        "phi": flat("SimTrack_boundaryMomentum_phi")[keep].astype(float),
        "pdg": flat("SimTrack_pdgId")[keep],
        "ekin": np.nan,
    })


def caloparticle_frame(tree):
    """One row per CaloParticle: (event, cp, pdg, energy, eta, phi, track)."""
    cols = {b: tree[f"CaloPart_{b}"].array(library="np")
            for b in ("pdgId", "energy", "eta", "phi", "trackId")}
    counts = np.array([len(x) for x in cols["pdgId"]])
    flat = {b: np.concatenate(v) for b, v in cols.items()}
    return pd.DataFrame({
        "event": np.repeat(np.arange(len(counts)), counts),
        "cp": np.concatenate([np.arange(n) for n in counts]),
        "pdg": flat["pdgId"],
        "energy": flat["energy"].astype(float),
        "eta": flat["eta"].astype(float),
        "phi": flat["phi"].astype(float),
        "track": flat["trackId"],
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", type=Path, help="one or two nanoML files")
    ap.add_argument("--tag")
    args = ap.parse_args()
    if len(args.files) > 2:
        ap.error("one file (its partner is found) or two")

    for path in args.files:
        d = describe(path)
        print(f"{path.name}: {LABELS[d['config']]} ({d['config']}), {d['particle']} "
              f"{d['energy']}, {d['events']} events, tables {d['tables']}")
    s = load_sample(*args.files, tag=args.tag)
    print(f"\nsample {s.tag}: {s.particle} {s.energy}")
    for key in KEYS:
        print(f"  {LABELS[key]:9s} {s.files[key].name}")
    print(f"  table read by the standard run: {s.table}")
    print(f"  mapping: {'found' if existing_tag(*s.files.values()) else 'none yet'}")


if __name__ == "__main__":
    main()
