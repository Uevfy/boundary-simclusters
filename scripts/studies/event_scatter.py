"""One event's boundary SimClusters, drawn where they cross: the event display.

    python scripts/event_scatter.py bound_photon_fc0_200GeV_testNanoML.root --row 352
    python scripts/event_scatter.py a.root b.root --row 352 --calopart 0 --table legacy
    python scripts/studies/event_scatter.py 11_fc1_s1_testNanoML.root --row 102 --single

In the house style: one
panel per figure, open markers with a shape and colour per particle type, marker
size by energy, the species legend under the panel. Two views per CaloParticle and
config, with the same axis ranges in both configs of a pair, so the panels compare:

    scatter_xy   the transverse (x, y) footprint on the boundary face
    scatter_rz   the longitudinal (|z|, r) view, both configs' front faces marked

Given a Standard/FineCalo pair (one file and its partner is found), `--row` is
the Standard row and the FineCalo side is **the same physical event**: its row
and its CaloParticles come from the match step's mappings, since the files
are not row-aligned (docs/data_handling.md) and CaloParticle order differs in some
events. Marker sizes share one scale across both configs of the event, so
they can be compared. `--single` draws only the file given, at its own row,
without pairing.

Marker size is kinetic energy, not `boundaryEnergy`, which includes rest mass
(docs/data_handling.md). Everything drawn is read from the branches: no parent
links or production vertices exist (docs/data_handling.md), so nothing here is a
reconstructed shower history.

Figures in outputs/<particle>/<energy>/<table>/ as
scatter_{xy,rz}_<tag>_ev<row>_cp<n>_<config>.png, where <row> is the Standard
row (the file's own row with --single) and <n> the Standard CaloParticle index.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from boundary_simclusters import samples  # noqa: E402
from boundary_simclusters.boundary_geometry import kinetic_energy  # noqa: E402
from boundary_simclusters.hep_style import BASE_FONT, INK_2, save, use_house_style  # noqa: E402
from boundary_simclusters.output_paths import table_dir  # noqa: E402

LABELS = samples.LABELS
SURFACES = {"fc0": 293.5, "fc1": 321.05}      # |z| of each config's front face, cm
PARTICLE_TEX = {"electron": "$e^-$", "positron": "$e^+$", "photon": r"$\gamma$"}

# PDG ID to label, colour and marker shape. Colour here encodes the particle
# type, not the config: each figure shows one config, named in its title.
PDG_STYLE_MAP = {
    22: {"label": r"$\gamma$", "color": "#1f77b4", "marker": "o"},
    11: {"label": r"$e^-$", "color": "#ff7f0e", "marker": "s"},
    -11: {"label": r"$e^+$", "color": "#2ca02c", "marker": "*"},
    13: {"label": r"$\mu^-$", "color": "#8c564b", "marker": ">"},
    -13: {"label": r"$\mu^+$", "color": "#e377c2", "marker": "<"},
    211: {"label": r"$\pi^+$", "color": "#d62728", "marker": "^"},
    -211: {"label": r"$\pi^-$", "color": "#9467bd", "marker": "v"},
    2112: {"label": r"$n$", "color": "#17becf", "marker": "h"},
    2212: {"label": r"$p$", "color": "#bcbd22", "marker": "D"},
}
DEFAULT_STYLE = {"label": "Other", "color": "#7f7f7f", "marker": "p"}
# Plain names for printed output, where mathtext would show as markup.
PLAIN = {22: "photon", 11: "e-", -11: "e+", 13: "mu-", -13: "mu+", 211: "pi+", -211: "pi-",
         2112: "n", 2212: "p"}


def style(pdg):
    return PDG_STYLE_MAP.get(int(pdg), {"label": f"PDG {int(pdg)}", **{
        k: DEFAULT_STYLE[k] for k in ("color", "marker")}})


# --- data -------------------------------------------------------------------------


def load_event(tree, row, table):
    """One event's SimClusters of `table` and its CaloParticles, as arrays."""
    p = samples.TABLES[table]
    names = [f"{p}_{b}" for b in ("impactPoint_x", "impactPoint_y", "impactPoint_z",
                                  "boundaryEnergy", "boundaryPmag", "pdgId", "CaloPartIdx")]
    names += ["CaloPart_pdgId", "CaloPart_energy", "CaloPart_eta"]
    a = {k: np.asarray(v[0]) for k, v in
         tree.arrays(names, entry_start=row, entry_stop=row + 1, library="np").items()}
    x, y, z = (a[f"{p}_impactPoint_{c}"].astype(float) for c in "xyz")
    return {
        "x": x, "y": y, "z": z, "abs_z": np.abs(z), "r": np.hypot(x, y),
        "energy": kinetic_energy(a[f"{p}_boundaryEnergy"].astype(float),
                                 a[f"{p}_boundaryPmag"].astype(float)),
        "pdg": a[f"{p}_pdgId"],
        "calopart": a[f"{p}_CaloPartIdx"],
        "caloparts": {b: a[f"CaloPart_{b}"] for b in ("pdgId", "energy", "eta")},
    }


def select_calopart(hits_data, calopart):
    """Restricts hits to one CaloParticle, or unchanged if None."""
    if calopart is None:
        return hits_data
    mask = hits_data["calopart"] == calopart
    out = {k: v[mask] for k, v in hits_data.items() if k != "caloparts"}
    out["caloparts"] = hits_data["caloparts"]
    return out


def calopart_label(hits_data, calopart):
    """The CaloParticle part of a figure title."""
    cp = hits_data["caloparts"]
    return (f" — CaloParticle {calopart} ({style(cp['pdgId'][calopart])['label']}, "
            f"{cp['energy'][calopart]:.0f} GeV, $\\eta$ = {cp['eta'][calopart]:+.2f})")


# --- drawing ----------------------------------------------------


def compute_point_sizes(energy, min_size=35, max_size=320, ref=None):
    """Scales marker size linearly with energy, against the reference range."""
    ref = energy if ref is None else ref
    if len(energy) == 0 or len(ref) == 0 or np.max(ref) == np.min(ref):
        return np.full_like(energy, min_size)
    norm = (energy - np.min(ref)) / (np.max(ref) - np.min(ref) + 1e-6)
    return min_size + norm * (max_size - min_size)


def padded_limits(lo, hi, frac=0.15, min_pad=5.0):
    """Widened axis limits to prevent edge markers from clipping."""
    pad = max((hi - lo) * frac, min_pad)
    return lo - pad, hi + pad


def species_legend(fig, pdg, legend_title, ncol=8):
    """Legend placed below the panel with matching symbols."""
    handles, labels = [], []
    for p in np.unique(pdg):
        cfg = style(p)
        handles.append(Line2D([], [], marker=cfg["marker"], linestyle="none", markersize=9,
                              color=cfg["color"], markeredgecolor=cfg["color"],
                              markerfacecolor="none", markeredgewidth=1.8))
        labels.append(cfg["label"])
    fig.legend(handles, labels, title=legend_title, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), ncol=min(ncol, len(handles)), frameon=False)


def draw(hits_data, xkey, ykey, title, xlabel, ylabel, figsize, ref_energy, path, legend_title,
         limits, surfaces=False):
    """One view of one CaloParticle in one config, saved to `path`."""
    fig, ax = plt.subplots(figsize=figsize)
    sizes = compute_point_sizes(hits_data["energy"], ref=ref_energy)
    # Largest first, so smaller markers are drawn on top of them.
    for i in np.argsort(hits_data["energy"])[::-1]:
        cfg = style(hits_data["pdg"][i])
        ax.scatter(hits_data[xkey][i], hits_data[ykey][i], s=sizes[i], facecolors="none",
                   edgecolors=cfg["color"], marker=cfg["marker"], linewidths=1.8, alpha=0.85)
    ax.set_title(title, pad=12, fontsize=BASE_FONT - 10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    if surfaces:
        for key, z in SURFACES.items():
            ax.axvline(z, ymax=0.86, color=INK_2, linestyle=":", linewidth=1.4)   # stops below the label
            ax.annotate(f"{LABELS[key]}'s surface", xy=(z, 0.03), xycoords=("data", "axes fraction"),
                        xytext=(-5, 0), textcoords="offset points", rotation=90, ha="right",
                        va="bottom", fontsize=BASE_FONT - 7, color=INK_2)
    species_legend(fig, hits_data["pdg"], legend_title)
    fig.tight_layout(rect=(0.0, 0.14, 1.0, 0.96))
    # Widen to fit the long title.
    return save(fig, path)


def view_limits(sels, xkey, ykey, surfaces, rmin0):
    """The axis ranges of one view, the same for every config of the CaloParticle."""
    xs = np.concatenate([s[xkey] for s in sels])
    ys = np.concatenate([s[ykey] for s in sels])
    if surfaces:                                  # the faces stay in range
        xs = np.concatenate([xs, list(SURFACES.values())])
    lo, hi = padded_limits(ys.min(), ys.max())
    return padded_limits(xs.min(), xs.max()), (max(0.0, lo) if rmin0 else lo, hi)


def draw_views(hits, calopart, head, ref_energy, stem, who, legend_title, others):
    """Both views of one CaloParticle; skipped, with a note, if it has no SimCluster."""
    sel = select_calopart(hits, calopart)
    title = f"{head}{calopart_label(hits, calopart)}"
    counts = ", ".join(f"{n} {PLAIN.get(int(p), int(p))}"
                       for p, n in zip(*np.unique(sel["pdg"], return_counts=True)))
    print(f"  {who}, CaloParticle {calopart}: "
          f"{len(sel['x'])} SimClusters{f' ({counts})' if counts else ''}")
    if not len(sel["x"]):
        print("    no boundary SimClusters, no figures")
        return
    for view, (xk, yk, name, xl, yl, size, rz) in {
        "xy": ("x", "y", "Front Face $(x, y)$", r"$x$ of impact point [cm]",
               r"$y$ of impact point [cm]", (7, 7.5), False),
        "rz": ("abs_z", "r", "Longitudinal $(|z|, r)$", r"$|z|$ of impact point [cm]",
               r"$r$ of impact point [cm]", (7.5, 7.5), True),
    }.items():
        path = stem.with_name(stem.name.format(view=view))
        limits = view_limits([sel] + others, xk, yk, rz, rz)
        print("    wrote", draw(sel, xk, yk, f"{title} — {name}", xl, yl, size,
                               ref_energy, path, legend_title, limits, rz))


# --- driver ------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", type=Path,
                    help="one nanoML file (its partner is found) or two")
    ap.add_argument("--row", type=int, required=True,
                    help="the event: the Standard row of a pair, or the file's own row with --single")
    ap.add_argument("--calopart", type=int,
                    help="one CaloParticle (Standard index); default: each that has SimClusters")
    ap.add_argument("--table", choices=samples.TABLES,
                    help="default: boundary where the file(s) carry it, else legacy")
    ap.add_argument("--tag", help="filename tag; default: the sample's (the file stem with --single)")
    ap.add_argument("--single", action="store_true",
                    help="draw only the one file given, at its own row, without pairing")
    args = ap.parse_args()
    if len(args.files) > 2 or (args.single and len(args.files) != 1):
        ap.error("one or two files; --single takes exactly one")

    if args.single:
        d = samples.describe(args.files[0])
        table = args.table or d["tables"][0]
        tree = samples.events(args.files[0])
        rows, cps = {d["config"]: args.row}, {d["config"]: None}
        trees, particle, energy = {d["config"]: tree}, d["particle"], d["energy"]
        tag = args.tag or args.files[0].stem
    else:
        sample = samples.load_sample(*args.files, tag=args.tag)
        table = args.table or sample.table
        mapping, caloparts = sample.mappings()
        partner = mapping.loc[mapping.fc0_row == args.row, "fc1_row"]
        if partner.empty:
            raise SystemExit(f"Standard row {args.row} has no FineCalo partner in {sample.tag}")
        rows = {"fc0": args.row, "fc1": int(partner.iloc[0])}
        cp_map = caloparts[caloparts.fc0_row == args.row]
        cps = {"fc0": None, "fc1": dict(zip(cp_map.fc0_cp.astype(int), cp_map.fc1_cp.astype(int)))}
        trees = {k: sample.tree(k) for k in samples.KEYS}
        particle, energy, tag = sample.particle, sample.energy, sample.tag
        tree = trees["fc0"]
    if table not in samples.TABLES:
        raise SystemExit(f"unknown table {table}")
    tag += "_bound" if table == "boundary" else ""
    figdir = table_dir(tree, table)

    hits = {k: load_event(trees[k], rows[k], table) for k in trees}
    # One size scale across both configs of the event, so markers compare.
    ref_energy = np.concatenate([h["energy"] for h in hits.values()])
    # Each CaloParticle the matching pairs (or, with --single, each with SimClusters).
    if args.calopart is not None:
        wanted = [args.calopart]
    elif args.single:
        wanted = sorted(np.unique(next(iter(hits.values()))["calopart"]).tolist())
    else:
        wanted = sorted(cps["fc1"])
    sample_tex = f"{PARTICLE_TEX.get(particle, particle)} {energy.replace('GeV', ' GeV')}"
    print(f"{particle} {energy}, {table} table, rows "
          + ", ".join(f"{LABELS[k]} {rows[k]}" for k in hits))

    use_house_style()
    for cp in wanted:
        own_cp = {k: cp if cps[k] is None else cps[k][cp] for k in hits}
        sels = {k: select_calopart(hits[k], own_cp[k]) for k in hits}
        for key in hits:
            own = own_cp[key]
            head = f"{LABELS[key]}, {table}, {sample_tex} (event {rows[key]})"
            stem = figdir / f"scatter_{{view}}_{tag}_ev{args.row}_cp{cp}_{key}.png"
            draw_views(hits[key], own, head, ref_energy, stem, f"{LABELS[key]} (row {rows[key]})",
                       f"Particle type ({table} collection, {energy.replace('GeV', ' GeV')})",
                       [s for s in sels.values() if len(s["x"])])


if __name__ == "__main__":
    main()
