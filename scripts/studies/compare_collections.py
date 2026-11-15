"""The boundary SimCluster collection against the legacy one, on the same events.

    python scripts/studies/compare_collections.py bound_*_testNanoML.root

Takes the nanoML files of several samples, pairs them among themselves
(`samples.py`), and reads both SimCluster tables of every pair that carries
the boundary one: `BoundarySimCluster_*`
(`mix:MergedCaloTruthBoundaryTrackSimCluster`) and the legacy `SimCluster_*`
(`mix:MergedCaloTruth`). Both describe the same events in the same file, so no
event pairing is needed between them; the Standard/FineCalo CaloParticle pairing
is the match step's (`calopart_matching_<tag>.parquet`, never by index,
docs/data_handling.md). Run scripts/analyze.py on each pair first.

Per sample and per (config, collection):
- **SimClusters per CaloParticle** (total, mean, median), and how many legacy
  SimClusters one boundary SimCluster stands for
- **energy closure**: sum of kinetic energies per CaloParticle over the gun
  energy (kinetic, not `boundaryEnergy`, which includes rest mass)
- **species mix**, **faces** (main face, other faces, interior, bore entries) and the share of **soft** SimClusters (below 0.1% of the gun energy)
- **the Standard/FineCalo comparison on each collection**: ratio of mean
  counts, share of CaloParticles with equal counts, FineCalo more or fewer

Figures, one colour per (config, collection) -- purple Standard and gold
FineCalo on the boundary collection, blue Standard and red FineCalo on the
legacy one (`hep_style.CONFIG_COLLECTION`); the line style is the energy (solid
200 GeV, dashed 20 GeV), as in every figure: per sample in
outputs/<particle>/<energy>/collections/ (SimClusters per CaloParticle, kinetic
energy spectrum, energy closure); across samples in outputs/slides/collections/
(FineCalo/Standard ratio, agreeing CaloParticles, closure), the collection
encoded by filled (boundary) or open (legacy) markers.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from boundary_simclusters import samples
from boundary_simclusters.boundary_geometry import label_crossings
from boundary_simclusters.hep_style import (BASE_FONT, CONFIG_COLLECTION, INK, energy_linestyle,
                                            energy_text, legend_with_headroom, save, series_label,
                                            single_figure,
                                            use_house_style)
from boundary_simclusters.output_paths import slide_path, style_dir

KEYS, LABELS = samples.KEYS, samples.LABELS
COLLECTIONS = ("boundary", "legacy")
FILL = {"boundary": "full", "legacy": "none"}
PARTICLE_TEX = {"electron": "$e^-$", "photon": r"$\gamma$"}
SOFT = 1e-3                                    # of the gun energy
HADRONS = (2112, 2212, 211, 321, 130, 310, 3122)


# --- reading ---------------------------------------------------------------------


def crossings(sample, table, key):
    d = label_crossings(samples.simcluster_frame(sample.tree(key), table, key), key)
    bore = (d["group"] == "inner side") & d["entering"] if key == "fc1" else np.zeros(len(d), bool)
    return d.assign(bore=bore)


def per_caloparticle(frames, caloparts):
    """{(collection, key): DataFrame with n and sum of E_kin per matched CaloParticle}."""
    rows = {"fc0": list(zip(caloparts.fc0_row, caloparts.fc0_cp)),
            "fc1": list(zip(caloparts.fc1_row, caloparts.fc1_cp))}
    out = {}
    for (collection, key), d in frames.items():
        g = d.groupby(["event", "cp"]).agg(n=("ekin", "size"), ekin=("ekin", "sum"))
        out[(collection, key)] = g.reindex(rows[key]).fillna(0).reset_index(drop=True)
    return out


def load(path_files):
    fc0s = [f for f in path_files if samples.describe(f)["config"] == "fc0"]
    data = {}
    for f in fc0s:
        sample = samples.load_sample(f, samples.find_partner(f, candidates=path_files))
        if "boundary" not in sample.tables:
            print(f"{sample.tag}: no boundary table, skipped")
            continue
        _, caloparts = sample.mappings()
        gun = float(np.median(samples.caloparticle_frame(sample.tree("fc0")).energy))
        frames = {(c, k): crossings(sample, c, k) for c in COLLECTIONS for k in KEYS}
        data[sample.tag] = dict(sample=sample, gun=gun, frames=frames,
                                cp=per_caloparticle(frames, caloparts))
    order = lambda t: (data[t]["sample"].particle, -data[t]["gun"])
    return {t: data[t] for t in sorted(data, key=order)}


# --- numbers ---------------------------------------------------------------------


def describe(s):
    """Rows per (config, collection) for one sample."""
    rows = {}
    for key in KEYS:
        for c in COLLECTIONS:
            d, cp = s["frames"][(c, key)], s["cp"][(c, key)]
            pdg = np.abs(d.pdg.to_numpy())
            closure = cp.ekin / s["gun"]
            rows[f"{LABELS[key]}, {c}"] = {
                "SimClusters (file)": len(d),
                "per CP mean": round(cp.n.mean(), 2),
                "per CP median": cp.n.median(),
                "sum Ekin/E_gun mean": round(closure.mean(), 3),
                "sum Ekin/E_gun median": round(closure.median(), 3),
                "photons %": round(100 * (pdg == 22).mean(), 1),
                "e+- %": round(100 * (pdg == 11).mean(), 1),
                "hadrons %": round(100 * np.isin(pdg, HADRONS).mean(), 1),
                "main face %": round(100 * (d["class"] == "main face").mean(), 1),
                "interior": int((d["class"] == "interior").sum()),
                "bore entries": int(d.bore.sum()),
                f"soft (<{SOFT:g} E_gun) %": round(100 * (d.ekin < SOFT * s["gun"]).mean(), 1),
            }
    return pd.DataFrame.from_dict(rows, orient="index")   # one row per (config, collection)


def conclusions(s):
    """The Standard/FineCalo comparison on each collection."""
    out = {}
    for c in COLLECTIONS:
        n0, n1 = s["cp"][(c, "fc0")].n, s["cp"][(c, "fc1")].n
        out[c] = {"ratio of means": round(n1.mean() / n0.mean(), 2),
                  "equal count %": round(100 * (n0 == n1).mean(), 1),
                  "FineCalo more %": round(100 * (n1 > n0).mean(), 1),
                  "FineCalo fewer %": round(100 * (n1 < n0).mean(), 1),
                  "closure Std": round((s["cp"][(c, "fc0")].ekin / s["gun"]).mean(), 3),
                  "closure Fine": round((s["cp"][(c, "fc1")].ekin / s["gun"]).mean(), 3)}
    return pd.DataFrame(out)


def legacy_per_boundary(s):
    """How many legacy SimClusters one boundary SimCluster stands for, per config."""
    return {LABELS[k]: round(len(s["frames"][("legacy", k)]) / len(s["frames"][("boundary", k)]), 1)
            for k in KEYS}


def label(s):
    return f"{s['sample'].particle} {s['sample'].energy.replace('GeV', ' GeV')}"


def tex(s):
    return f"{PARTICLE_TEX.get(s['sample'].particle, s['sample'].particle)} " \
           f"{s['sample'].energy.replace('GeV', ' GeV')}"


# --- per-sample figures ----------------------------------------------------------


def series(ax, values, bins, key, collection, dash, energy, norm=1.0, per_width=False):
    counts, _ = np.histogram(values, bins=bins)
    counts = counts / norm / (np.diff(bins) if per_width else 1.0)
    ax.stairs(counts, bins, color=CONFIG_COLLECTION[(key, collection)], linestyle=dash,
              linewidth=2, label=series_label(key, collection, energy))


def finish(fig, ax, path, title, xlabel, ylabel, loc="upper right"):
    ax.set_xlabel(xlabel, linespacing=1.5)
    ax.set_ylabel(ylabel, linespacing=1.5)
    ax.set_title(title, pad=12)
    legend_with_headroom(ax, loc=loc, fontsize=BASE_FONT - 3)
    return save(fig, path)


def figures_sample(s):
    folder = style_dir(s["sample"].tree("fc0"), "collections")
    tag = s["sample"].tag
    energy = s["sample"].energy
    dash = energy_linestyle(energy)                    # solid 200 GeV, dashed 20 GeV
    written = []

    fig, ax = single_figure()
    # Log axis: integer bins up to 10, then wider ones, so 1 and 150 both read.
    # CaloParticles with no SimCluster (a few per sample) fall off the axis.
    top = max(int(cp.n.max()) for cp in s["cp"].values())
    bins = np.unique(np.r_[np.arange(0.5, 10.5), np.round(np.geomspace(10.5, top + 1.5, 25)) + 0.5])
    for key in KEYS:
        for c in COLLECTIONS:
            series(ax, s["cp"][(c, key)].n, bins, key, c, dash, energy, per_width=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, bins[-1])
    written.append(finish(fig, ax, folder / f"n_per_caloparticle_{tag}.png",
                          f"SimClusters per CaloParticle, {tex(s)}",
                          "SimClusters per CaloParticle",
                          "CaloParticles per unit of count"))

    fig, ax = single_figure()
    bins = np.logspace(-6, np.log10(1.1), 61)
    n_cp = len(s["cp"][("boundary", "fc0")])
    for key in KEYS:
        for c in COLLECTIONS:
            series(ax, s["frames"][(c, key)].ekin / s["gun"], bins, key, c, dash, energy, norm=n_cp)
    ax.set_xscale("log")
    ax.set_yscale("log")
    written.append(finish(fig, ax, folder / f"ekin_spectrum_{tag}.png",
                          f"Kinetic energy of each SimCluster, {tex(s)}",
                          r"$E_\mathrm{kin}$ / $E_\mathrm{gun}$", "SimClusters per CaloParticle per bin",
                          loc="upper left"))

    fig, ax = single_figure()
    bins = np.linspace(0, 1.1, 56)
    for key in KEYS:
        for c in COLLECTIONS:
            series(ax, s["cp"][(c, key)].ekin / s["gun"], bins, key, c, dash, energy)
    ax.set_yscale("log")
    written.append(finish(fig, ax, folder / f"closure_{tag}.png",
                          f"Kinetic energy recorded per CaloParticle, {tex(s)}",
                          r"$\Sigma E_\mathrm{kin}$ / $E_\mathrm{gun}$ per CaloParticle", "CaloParticles",
                          loc="upper left"))
    return written


# --- cross-sample slides ---------------------------------------------------------


def slide(data, values, name, title, ylabel, per_config=False, fmt="{:.2f}"):
    """`values(s, collection[, key])` per sample; filled markers boundary, open legacy."""
    use_house_style()
    fig, ax = single_figure()
    tags = list(data)
    x = np.arange(len(tags))
    groups = [(c, k) for c in COLLECTIONS for k in KEYS] if per_config else [(c, None) for c in COLLECTIONS]
    offsets = np.linspace(-0.18, 0.18, len(groups))
    for (c, k), dx in zip(groups, offsets):
        y = [values(data[t], c, k) if per_config else values(data[t], c) for t in tags]
        colour = CONFIG_COLLECTION[(k, c)] if per_config else INK
        energies = [data[t]["sample"].energy for t in tags]
        name_ = series_label(k, c, energies) if per_config else f"{c}, {energy_text(energies)}"
        ax.plot(x + dx, y, linestyle="none", marker="o", markersize=10, color=colour,
                fillstyle=FILL[c], markeredgewidth=2, label=name_)
        if not per_config:
            for xi, yi in zip(x + dx, y):
                ax.annotate(fmt.format(yi), (xi, yi), textcoords="offset points", xytext=(12, -4),
                            fontsize=BASE_FONT - 4, color=INK)
    ax.set_xticks(x, [tex(data[t]) for t in tags])
    ax.set_xlim(-0.6, len(tags) - 0.4)
    ax.set_ylim(bottom=0)
    ax.set_ylabel(ylabel, linespacing=1.5)
    ax.set_title(title, pad=12)
    legend_with_headroom(ax, loc="upper right", fontsize=BASE_FONT - 3)
    return save(fig, slide_path(name, "collections"))


def slides(data):
    def ratio(s, c):
        return s["cp"][(c, "fc1")].n.mean() / s["cp"][(c, "fc0")].n.mean()

    def equal(s, c):
        return 100 * (s["cp"][(c, "fc0")].n == s["cp"][(c, "fc1")].n).mean()

    def closure(s, c, k):
        return (s["cp"][(c, k)].ekin / s["gun"]).mean()

    return [slide(data, ratio, "ratio_by_collection.png",
                  "FineCalo / Standard SimClusters per CaloParticle", "ratio of means"),
            slide(data, equal, "agreement_by_collection.png",
                  "CaloParticles with the same count in both configs", "CaloParticles  [%]",
                  fmt="{:.1f}%"),
            slide(data, closure, "closure_by_collection.png",
                  "Kinetic energy recorded per CaloParticle",
                  r"mean $\Sigma E_\mathrm{kin}$ / $E_\mathrm{gun}$", per_config=True)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", type=Path, help="nanoML files of samples with the boundary table")
    args = ap.parse_args()
    data = load(args.files)
    if not data:
        raise SystemExit("no sample with a boundary table among the files")
    use_house_style()
    pd.set_option("display.width", 200)
    for tag, s in data.items():
        print(f"\n=== {label(s)} ({tag}), {len(s['cp'][('boundary', 'fc0')])} matched CaloParticles ===")
        print(describe(s).to_string())
        print("\nStandard vs FineCalo on each collection:")
        print(conclusions(s).to_string())
        print("legacy SimClusters per boundary SimCluster:", legacy_per_boundary(s))
        for path in figures_sample(s):
            print("wrote", path)
    print()
    for path in slides(data):
        print("wrote", path)


if __name__ == "__main__":
    main()
