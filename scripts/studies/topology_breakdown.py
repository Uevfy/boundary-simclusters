"""The decomposition classes broken down by species, energy and |eta|, across samples.

    python scripts/studies/topology_breakdown.py bound_*_testNanoML.root

Takes the nanoML files of several samples, pairs them among themselves
(`samples.py`) and reads each pair's classification
(`decompositions_<tag>_bound.parquet`, `decomposition_simclusters_<tag>_bound.parquet`,
written by the decompositions step of scripts/analyze.py). Species
and energy are the samples themselves; within a sample the breakdown is in the
CaloParticle's |eta|, in five bins over the gun band. Shares carry binomial
errors, sqrt(p(1 - p) / n).

Prints, per sample: the classes, and the classes per
|eta| bin -- beside the share of CaloParticles whose Standard side has a
charged SimCluster arriving below 1 GeV, which is what drives "other replacement"
-- and the topologies expected to concentrate
disagreements, each mapped onto what the boundary table records. Draws two
slides into outputs/slides/boundary/: lost in the gap vs |eta| (every sample)
and gap conversions vs |eta| (photon samples). Colour is the species (green electrons,
magenta photons, `hep_style.ELECTRON_COLOUR` / `PHOTON_COLOUR`: hues no config
scheme uses), line style the energy.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from boundary_simclusters import output_paths, samples
from boundary_simclusters.hep_style import (BASE_FONT, ELECTRON_COLOUR, PHOTON_COLOUR,
                                            energy_linestyle, legend_with_headroom, save, single_figure,
                                            use_house_style)
from boundary_simclusters.output_paths import slide_path
from boundary_simclusters.pipeline.classify_decompositions import CLASSES

ETA_BINS = np.linspace(1.70, 2.70, 6)          # the gun band, from GenPart_eta
SPECIES = {"electron": (ELECTRON_COLOUR, "$e^-$"), "photon": (PHOTON_COLOUR, r"$\gamma$")}
MARKER = {"200GeV": "o", "20GeV": "s"}
HADRONS = (2112, 2212, 211, 321, 130, 310, 3122)
CHARGED = (11, 13, 211, 321, 2212)
SLOW = 1.0      # GeV of kinetic energy at Standard's surface: "arrives nearly stopped"


def load(files):
    """{tag: (sample, classification, matching)} for every pair among `files`."""
    fc0s = [f for f in files if samples.describe(f)["config"] == "fc0"]
    out = {}
    for f in fc0s:
        sample = samples.load_sample(f, samples.find_partner(f, candidates=files))
        if "boundary" not in sample.tables:
            print(f"{sample.tag}: no boundary table, skipped")
            continue
        paths = [output_paths.PROCESSED / f"{name}_{sample.tag}_bound.parquet"
                 for name in ("decompositions", "decomposition_simclusters")]
        if not all(p.exists() for p in paths):
            raise SystemExit(f"{sample.tag}: run python scripts/analyze.py {Path(f).name} first "
                             "(needs the decompositions step)")
        out[sample.tag] = (sample, *(pd.read_parquet(p) for p in paths))
    order = lambda tag: (out[tag][0].particle, -float(out[tag][0].energy.replace("GeV", "")))
    return {tag: out[tag] for tag in sorted(out, key=order)}


def share(mask):
    """(percent, binomial error in percent) of a boolean array."""
    n = len(mask)
    p = mask.mean() if n else np.nan
    return 100 * p, 100 * np.sqrt(p * (1 - p) / n) if n else np.nan


def label(sample):
    return f"{sample.particle} {sample.energy.replace('GeV', ' GeV')}"


def report_classes(data):
    print("\nClasses by sample (share of paired CaloParticles)")
    print("-------------------------------------------------")
    rows = {}
    for tag, (sample, cls, _) in data.items():
        rows[label(sample)] = {c: f"{share(cls['class'] == c)[0]:6.2f}%" for c in CLASSES} | {"n": len(cls)}
    print(pd.DataFrame(rows).to_string())


def slow_arrivals(cls, match):
    """Per CaloParticle: does Standard record a charged SimCluster below `SLOW` GeV?"""
    std = match[match.config == "fc0"]
    slow = std[np.isin(np.abs(std.pdg), CHARGED) & (std.ekin < SLOW)].calopart.unique()
    return cls.calopart.isin(slow)


def report_eta(data):
    for tag, (sample, cls, match) in data.items():
        cls = cls.assign(slow=slow_arrivals(cls, match).to_numpy())
        bins = pd.cut(cls.abs_eta, ETA_BINS)
        print(f"\n{label(sample)}: by |eta| of the CaloParticle (share +- binomial error)")
        print("-" * (len(label(sample)) + 52))
        rows = {}
        for interval, group in cls.groupby(bins, observed=True):
            row = {"n": len(group)}
            for c in CLASSES[:5]:
                p, e = share(group["class"] == c)
                row[c] = f"{p:5.1f} +- {e:3.1f}"
            p, e = share(group.slow.to_numpy())
            row[f"Std charged < {SLOW:g} GeV"] = f"{p:5.1f} +- {e:3.1f}"
            rows[f"{interval.left:.1f}-{interval.right:.1f}"] = row
        print(pd.DataFrame(rows).T.to_string())


def checks(sample, key):
    """What the SimClusters of one config say about re-entry, backscattering and non-EM species."""
    tree = sample.tree(key)
    boundary = samples.simcluster_frame(tree, "boundary", key)
    legacy = samples.simcluster_frame(tree, "legacy", key)
    upstream = boundary[boundary["where"] == "main face"]
    backwards = int((np.sign(upstream.eta) != np.sign(upstream.z)).sum())
    twice = int((legacy.groupby(["event", "root_track"]).size() > 1).sum())
    other = boundary[boundary["where"] == "other face"]
    side = int((~other.face.str.startswith("inner bore")).sum())
    non_em = int((~np.isin(np.abs(boundary.pdg), (11, 22))).sum())
    hadrons = {name: f[np.isin(np.abs(f.pdg), HADRONS)] for name, f in
               (("boundary", boundary), ("legacy", legacy))}
    counted = "; ".join(f"{name} hadrons {len(h)} in {h.groupby(['event', 'cp']).ngroups} CaloParticles"
                        for name, h in hadrons.items())
    return (f"tracks with two legacy SimClusters: {twice}; moving backwards on the upstream face: "
            f"{backwards} of {len(upstream)}; on other faces, not bore: {side}; "
            f"not e or photon: {non_em} of {len(boundary)}; {counted}")


def uniformity_in_phi(sample, cls):
    """Rate of the two main causes in six bins of the CaloParticle's phi, against a flat rate."""
    phi = sample.tree("fc0")["CaloPart_phi"].array(library="np")
    angle = np.array([phi[row][cp] for row, cp in zip(cls.fc0_row.astype(int), cls.fc0_cp.astype(int))])
    bins = pd.cut(angle, np.linspace(-np.pi, np.pi, 7))
    out = []
    for cause in ("other replacement", "photon conversion"):
        hit = cls["class"].to_numpy() == cause
        if hit.sum() < 5:
            continue
        g = pd.Series(hit).groupby(bins, observed=False).agg(["mean", "size"])
        error = np.sqrt(g["mean"] * (1 - g["mean"]) / g["size"])
        out.append(f"{cause} chi2 {(((g['mean'] - hit.mean()) / error) ** 2).sum():.1f}/5")
    return "; ".join(out)


def report_topologies(data):
    """The topologies expected to concentrate disagreements."""
    print("\nExpected topologies, as the boundary table records them")
    print("-------------------------------------------------------")
    for tag, (sample, cls, match) in data.items():
        n = len(cls)
        replaced = cls[cls["class"] == "other replacement"].calopart
        members = match[(match.config == "fc1") & match.calopart.isin(replaced) & (match.role == "descendant")]
        hadronic = members[np.isin(np.abs(members.pdg), HADRONS)].calopart.nunique()
        fine = match[match.config == "fc1"]
        bore = fine[(fine["group"] == "inner side") & fine.entering]
        print(f"{label(sample)}  ({n} CaloParticles)")
        print(f"  photon conversions near the boundary  {100 * (cls['class'] == 'photon conversion').mean():5.2f}%")
        print(f"  showers / hadronic interactions       other replacement {100 * (cls['class'] == 'other replacement').mean():5.2f}%, "
              f"of them with a hadron among FineCalo's SimClusters: {hadronic}")
        print(f"  low-energy secondaries                no descendant {100 * (cls['class'] == 'no descendant').mean():5.2f}%, "
              f"mixed {100 * (cls['class'] == 'mixed').mean():5.2f}%")
        print(f"  side entries                          {len(bore)} in {bore.calopart.nunique()} CaloParticles, "
              f"by class {cls[cls.n_side_entries > 0]['class'].value_counts().to_dict()}")
        print("  leaving and re-entering              cannot add a FineCalo SimCluster (a track is flagged once)")
        for key in samples.KEYS:
            print(f"  {samples.LABELS[key]:9s} {checks(sample, key)}")
        print(f"  phi, six bins: {uniformity_in_phi(sample, cls)}")


def figure_vs_eta(data, cause, particles, name, title, ylabel):
    """Share of CaloParticles of class `cause`, per |eta| bin, one line per sample."""
    use_house_style()
    fig, ax = single_figure()
    centres = 0.5 * (ETA_BINS[:-1] + ETA_BINS[1:])
    shown = [tag for tag, (s, _, _) in data.items() if s.particle in particles]
    # Shift each sample a little sideways, so overlapping error bars stay apart.
    dodge = dict(zip(shown, 0.012 * (np.arange(len(shown)) - (len(shown) - 1) / 2)))
    for tag in shown:
        sample, cls, _ = data[tag]
        colour, tex = SPECIES.get(sample.particle, ("#9c9ca1", sample.particle))
        dash, marker = energy_linestyle(sample.energy), MARKER.get(sample.energy, "^")
        bins = pd.cut(cls.abs_eta, ETA_BINS)
        values = [share(g["class"] == cause) for _, g in cls.groupby(bins, observed=False)]
        p, e = np.array(values).T
        total = share(cls["class"] == cause)[0]
        drawn = ax.errorbar(centres + dodge[tag], p, yerr=e, color=colour, linestyle=dash,
                            linewidth=2, marker=marker, markersize=8, capsize=3,
                            label=f"{tex} {sample.energy.replace('GeV', ' GeV')}, boundary  "
                                  f"({total:.1f}% overall)")
        for bar in drawn[2]:
            bar.set_linestyle("-")              # error bars solid, whatever the line style
    ax.set_xlim(ETA_BINS[0], ETA_BINS[-1])
    ax.set_ylim(bottom=0)
    ax.set_xlabel("|η| of the CaloParticle", linespacing=1.5)
    ax.set_ylabel(ylabel, linespacing=1.5)
    ax.set_title(title, pad=12)
    legend_with_headroom(ax, loc="upper left", fontsize=BASE_FONT - 3)
    return save(fig, slide_path(name, "boundary"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", type=Path, help="nanoML files of several samples")
    args = ap.parse_args()
    data = load(args.files)
    if not data:
        raise SystemExit("no sample with a boundary table among the files")
    report_classes(data)
    report_eta(data)
    report_topologies(data)
    print()
    print("wrote", figure_vs_eta(data, "other replacement", SPECIES, "lost_in_gap_vs_eta_bound.png",
                                 "Reaching Standard's surface but not FineCalo's",
                                 "CaloParticles with an object\nlost in the gap  [%]"))
    if any(s.particle == "photon" for s, _, _ in data.values()):
        print("wrote", figure_vs_eta(data, "photon conversion", ("photon",),
                                     "gap_conversions_vs_eta_bound.png",
                                     "Photons converting between the two surfaces",
                                     "CaloParticles with a gap conversion  [%]"))


if __name__ == "__main__":
    main()
