"""Which face of its config's boundary volume each crossing hit: file in, census out.

    python scripts/analyze.py bound_photon_fc0_200GeV_testNanoML.root --steps surfaces
    python scripts/analyze.py a.root b.root --steps surfaces --table legacy --excel

The surfaces step of scripts/analyze.py, `run(sample)`.

The two mechanisms test two different Geant4 volumes -- Standard the `CALO`
polycone, FineCalo the `HGCal` polycone -- so a crossing is
classified by the face of its own config's volume it landed on, using the
volumes as CMSSW defines them (`boundary_geometry.classify_faces`). CMSSW itself
names only the volumes; the face names are ours, read off the polycone planes
in order.

Every crossing falls in one **class**:

    main face      the face nearly every forward track crosses: CALO's
                   tracker-cavity step (|z| = 293.5 cm), HGCal's front face
                   (321.05 cm)
    other faces    any other face of the volume
    interior       inside the volume, on no face: a track flagged at its
                   production vertex (Standard only)
    exterior       on no face and outside the volume

and every face to one **group**, from its geometry rather than its name:

    main face        as above
    inner side       faces along z whose outward normal points at the beam:
                     the bores, and CALO's tracker-cavity wall
    outer side       faces along z whose normal points away: walls and cones
    rear and steps   the other faces across z: rear face, caps, bore steps

A crossing on a side face is **entering** when its momentum points into the
volume. Entering inner-side crossings are the bore entries,
which have no counterpart in the other config by geometry, not by mechanism.

Two collections are reported. SimTracks are every crossing -- the right sample
for geometry. SimClusters (from `--table`: the boundary table where both files
carry it) are the crossings that also left a hit, what the analysis
compares; the figures and the per-SimCluster output are of those.

Figures, in outputs/<particle>/<energy>/<table>/ (`output_paths.table_dir`):

    surface_classes_<tag>        SimClusters per class, both configs:
                                 upstream face, other surfaces,
                                 not on boundary volume surface (interior +
                                 exterior)
    surface_faces_<group>_<tag>  SimClusters per face, one figure per group
                                 that has any (inner side, outer side, rear and
                                 steps; the main face is one bar)
    surface_ke_<group>_<tag>     kinetic energy of the group's SimClusters

Data: data/processed/surface_<tag>.parquet, one row per SimCluster with its face,
group, class and entering flag -- what the classification step reads to separate
bore entries and interior objects before matching.
"""


import numpy as np
import pandas as pd

from boundary_simclusters import samples  # noqa: E402
from boundary_simclusters.boundary_geometry import (MAIN_FACE, SOURCE_ENVELOPES,  # noqa: E402
                                                   face_groups, label_crossings)
from boundary_simclusters.hep_style import (BASE_FONT, CONFIG_COLLECTION, GUIDE, INK_2,  # noqa: E402
                                            config_colours,
                                            energy_linestyle, series_label,
                       legend_with_headroom, save, single_figure, use_house_style)
from boundary_simclusters import output_paths
from boundary_simclusters.output_paths import style_dir, table_dir  # noqa: E402

# FineCalo first, as in every legend of the repo; colour is bound to config and table.
KEYS = ("fc1", "fc0")
# Comparison figures overlay both tables: keys (config, table), see `compare`.
COMPARE_KEYS = (("fc1", "boundary"), ("fc0", "boundary"), ("fc1", "legacy"), ("fc0", "legacy"))
LABELS = samples.LABELS
COLORS = config_colours("boundary")       # reset per run to the table drawn
SERIES = {}                               # legend names, set per run
LINESTYLE = "-"                           # the sample's energy, per run
CLASSES = ("main face", "other faces", "interior", "exterior")
# The bars of surface_classes. CLASSES
# stay the data labels (parquet); the figure only renames them and
# merges the two classes that lie on no face.
FIGURE_CLASSES = {"upstream face": ("main face",),
                  "other surfaces": ("other faces",),
                  "not on boundary\nvolume surface": ("interior", "exterior")}
GROUPS = ("main face", "inner side", "outer side", "rear and steps")
GROUP_TITLES = {"main face": "the main face", "inner side": "the inner side (bores, cavity wall)",
                "outer side": "the outer side (walls, cones)",
                "rear and steps": "the rear face, caps and steps"}
PARTICLE = {22: "photon", 11: "e-", -11: "e+", 2112: "neutron", 2212: "proton",
            211: "pi+", -211: "pi-", 13: "mu-", -13: "mu+", 130: "K0L", 321: "K+", -321: "K-"}


# --- classification --------------------------------------------------------------


# --- tables ---------------------------------------------------------------------------


def summary(frames):
    """One row per config: how many crossings fell in each class."""
    rows = []
    for key in samples.KEYS:
        f = frames[key]
        row = {"config": LABELS[key], "volume": SOURCE_ENVELOPES[key]["volume"],
               "crossings": len(f)}
        for c in CLASSES:
            n = int((f["class"] == c).sum())
            row[c], row[f"{c} %"] = n, round(100 * n / len(f), 3)
        row["bore entries"] = int(((f.group == "inner side") & f.entering).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def species(group, top=2):
    """The leading particle types of a set of crossings, e.g. 'photon 60%, e- 30%'."""
    counts = group.pdg.map(lambda p: PARTICLE.get(int(p), str(int(p)))).value_counts()
    return ", ".join(f"{name} {100 * n / len(group):.0f}%" for name, n in counts.head(top).items())


def by_face(frames):
    """One row per (config, face) holding any crossing."""
    rows = []
    for key in samples.KEYS:
        f = frames[key]
        for face, g in f.groupby("face", sort=False):
            rows.append({
                "config": LABELS[key], "face": face, "group": g.group.iloc[0],
                "n": len(g), "share %": round(100 * len(g) / len(f), 3),
                "entering %": round(100 * g.entering.mean(), 1),
                "median KE": (round(float(np.nanmedian(g.ekin)), 4)
                              if np.isfinite(g.ekin).any() else np.nan),
                "species": species(g),
            })
    return (pd.DataFrame(rows).sort_values(["config", "n"], ascending=[True, False])
            .reset_index(drop=True))


def write_workbook(path, sheets):
    """{sheet: table} as a fresh workbook (never appended to: a stale sheet would
    look current). Counts and percentages are separate numeric columns, so the
    sheet sorts and sums."""
    with pd.ExcelWriter(path, engine="openpyxl", mode="w") as xl:
        for name, table in sheets.items():
            table.to_excel(xl, sheet_name=name[:31], index=False)
            sheet = xl.sheets[name[:31]]
            for i, col in enumerate(table.columns, start=1):
                width = max(len(str(col)), *(len(str(v)) for v in table[col]))
                sheet.column_dimensions[sheet.cell(1, i).column_letter].width = width + 3
    return path


# --- figures ----------------------------------------------------------------------------


def finish(fig, ax, path, *, title, xlabel, ylabel, loc="upper right"):
    ax.set_xlabel(xlabel, linespacing=1.5)
    ax.set_ylabel(ylabel, linespacing=1.5)
    ax.set_title(title, pad=12)
    legend_with_headroom(ax, loc=loc)
    return save(fig, path)


def config_of(key):
    """The config of a series key: "fc1", or the first part of ("fc1", "legacy")."""
    return key if isinstance(key, str) else key[0]


def figure_classes(frames, figdir, tag, prefix):
    """SimClusters per FIGURE_CLASSES bar, both configs, log scale: the upstream
    face holds >99%."""
    fig, ax = single_figure()
    x = np.arange(len(FIGURE_CLASSES))
    width = 0.76 / len(KEYS)
    for i, key in enumerate(KEYS):
        counts = [int(frames[key]["class"].isin(classes).sum())
                  for classes in FIGURE_CLASSES.values()]
        offset = (i - (len(KEYS) - 1) / 2) * width
        ax.bar(x + offset, np.maximum(counts, 0.0), width * 0.94, color=COLORS[key],
               label=f"{SERIES[key]}  ({SOURCE_ENVELOPES[config_of(key)]['volume']}, "
                     f"{sum(counts):,} SimClusters)")
        for xi, n in zip(x, counts):
            ax.annotate(f"{n:,}", xy=(xi + offset, max(n, 0.6)), xytext=(0, 3),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=BASE_FONT - (4 if len(KEYS) <= 2 else 6), color=INK_2 if n else GUIDE)
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(list(FIGURE_CLASSES))
    return finish(fig, ax, figdir / f"surface_classes_{tag}.png",
                  title="Where the SimClusters sit",
                  xlabel="Class of the crossing, on each configuration's own volume",
                  ylabel="SimClusters")


def figure_faces(frames, group, figdir, tag):
    """SimClusters on each face of one group, each volume's faces in its own
    block: the volumes have different faces, so no face is paired across them."""
    rows = []
    for key in KEYS:
        counts = frames[key].face.value_counts()
        for name, g in face_groups(config_of(key)).items():
            if g == group:
                rows.append((name, key, int(counts.get(name, 0))))
    fig, ax = single_figure()
    y = np.arange(len(rows))
    for i, (name, key, n) in enumerate(rows):
        ax.barh(i, n, 0.62, color=COLORS[key])
        ax.annotate(f"{n:,}", xy=(n, i), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=BASE_FONT - 4, color=INK_2 if n else GUIDE)
    for key in KEYS:                                   # legend entries, one per volume
        total = sum(n for _, k, n in rows if k == key)
        ax.barh(np.nan, 0, color=COLORS[key],
                label=f"{SERIES[key]}  ({SOURCE_ENVELOPES[config_of(key)]['volume']}, {total:,})")
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlim(0, max(max(n for *_, n in rows), 1) * 6)
    ax.set_yticks(y)
    ax.set_yticklabels([name for name, _, _ in rows], fontsize=BASE_FONT - 2)
    ax.invert_yaxis()
    ax.set_xlabel("SimClusters", linespacing=1.5)
    ax.set_ylabel("Face", linespacing=1.5)
    ax.set_title(f"SimClusters on {GROUP_TITLES[group]}", pad=12)
    # A plain legend: the rows are fixed, so there is no axis range to clear.
    ax.legend(loc="lower right")
    return save(fig, figdir / f"surface_faces_{group.replace(' ', '_')}_{tag}.png")


def figure_ke(frames, group, e_gun, figdir, tag):
    """Kinetic energy of the SimClusters on one group of faces."""
    fig, ax = single_figure()
    bins = np.logspace(-5, np.log10(1.1 * e_gun), 61)
    for key in KEYS:
        f = frames[key]
        ke = f.ekin[(f.group == group).to_numpy()].to_numpy()
        ke = ke[ke > 0]
        counts, edges = np.histogram(ke, bins=bins)
        med = f", median {np.median(ke):.3g} GeV" if len(ke) else ""
        ax.stairs(counts, edges, color=COLORS[key], linestyle=LINESTYLE,
                  label=f"{SERIES[key]}  ({len(ke):,}{med})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.5)
    return finish(fig, ax, figdir / f"surface_ke_{group.replace(' ', '_')}_{tag}.png",
                  title=f"Kinetic energy of the SimClusters on {GROUP_TITLES[group]}",
                  xlabel="$E_\\mathrm{kin}$ at the crossing  [GeV]",
                  ylabel="SimClusters")


# --- driver ------------------------------------------------------------------------------


def run(sample, table=None, tag=None, excel=False):
    """The face census of `sample` on `table` (default: boundary where both files
    carry it): prints the tables, draws the figures, writes
    surface_<tag>.parquet, and returns the per-SimCluster labels it wrote."""
    global COLORS, SERIES, LINESTYLE
    table = table or sample.table
    if table not in sample.tables:
        raise SystemExit(f"{table} table not in both files ({sample.tables})")
    COLORS = config_colours(table)
    SERIES = {k: series_label(k, table, sample.energy) for k in samples.KEYS}
    LINESTYLE = energy_linestyle(sample.energy)
    tag = (tag or sample.tag) + ("_bound" if table == "boundary" else "")
    trees = {k: sample.tree(k) for k in samples.KEYS}
    figdir = table_dir(trees["fc0"], table)

    print(f"sample {sample.tag}: {sample.particle} {sample.energy}, {table} table")
    for key in samples.KEYS:
        print(f"  {LABELS[key]:9s} {sample.files[key].name}  "
              f"(volume {SOURCE_ENVELOPES[key]['volume']}, main face: {MAIN_FACE[key]})")

    clusters = {k: label_crossings(samples.simcluster_frame(trees[k], table, k), k) for k in samples.KEYS}
    tracks = {k: label_crossings(samples.simtracks(trees[k]), k) for k in samples.KEYS}

    sheets = {}
    for name, frames in (("SimTrack", tracks), ("SimCluster", clusters)):
        s, f = summary(frames), by_face(frames)
        print(f"\n{name}s by class\n{'-' * (len(name) + 10)}\n{s.to_string(index=False)}")
        print(f"\n{name}s by face\n{'-' * (len(name) + 9)}\n{f.to_string(index=False)}")
        sheets[f"{name} summary"], sheets[f"{name} by face"] = s, f

    use_house_style()
    prefix = f"{samples.TABLES[table]}_impactPoint"
    print("\nwrote", figure_classes(clusters, figdir, tag, prefix))
    e_gun = float(np.median(np.concatenate(trees["fc0"]["CaloPart_energy"].array(library="np"))))
    for group in GROUPS:
        n = sum(int((clusters[k].group == group).sum()) for k in samples.KEYS)
        if not n:
            print(f"no SimClusters on {GROUP_TITLES[group]}: no figures for it")
            continue
        if group != "main face":
            print("wrote", figure_faces(clusters, group, figdir, tag))
        print("wrote", figure_ke(clusters, group, e_gun, figdir, tag))

    out = output_paths.PROCESSED / f"surface_{tag}.parquet"
    cols = ["event", "cp", "pdg", "face", "group", "class", "entering", "ekin", "abs_z", "r"]
    labels = pd.concat([clusters[k][cols].assign(config=k) for k in samples.KEYS],
                       ignore_index=True)
    labels.to_parquet(out, index=False)
    print(f"wrote {out}")
    if excel:
        print("wrote", write_workbook(figdir / f"surfaces_{tag}.xlsx", sheets))
    return labels


def compare(sample, tag=None):
    """The surface figures with both SimCluster collections overlaid -- four
    series, one colour each -- into outputs/<particle>/<energy>/comparisons/
    <name>_<tag>_compare.png. Figures only; the census and the surface parquet
    stay per table (`run`). Returns the paths."""
    global KEYS, COLORS, SERIES, LINESTYLE
    if not {"boundary", "legacy"} <= set(sample.tables):
        raise SystemExit(f"{sample.tag} does not carry both SimCluster tables ({sample.tables})")
    saved = KEYS, COLORS, SERIES, LINESTYLE
    written = []
    try:
        KEYS = COMPARE_KEYS
        COLORS = {k: CONFIG_COLLECTION[k] for k in KEYS}
        SERIES = {k: series_label(k[0], k[1], sample.energy) for k in KEYS}
        LINESTYLE = energy_linestyle(sample.energy)
        trees = {k: sample.tree(k) for k in samples.KEYS}
        figdir = style_dir(trees["fc0"], "comparisons")
        tag = (tag or sample.tag) + "_compare"
        clusters = {k: label_crossings(samples.simcluster_frame(trees[k[0]], k[1], k[0]), k[0])
                    for k in KEYS}
        use_house_style()
        written.append(figure_classes(clusters, figdir, tag, "SimCluster_impactPoint"))
        e_gun = float(np.median(np.concatenate(trees["fc0"]["CaloPart_energy"].array(library="np"))))
        for group in GROUPS:
            if not sum(int((clusters[k].group == group).sum()) for k in KEYS):
                continue
            if group != "main face":
                written.append(figure_faces(clusters, group, figdir, tag))
            written.append(figure_ke(clusters, group, e_gun, figdir, tag))
    finally:
        KEYS, COLORS, SERIES, LINESTYLE = saved
    for path in written:
        print("wrote", path)
    return written
