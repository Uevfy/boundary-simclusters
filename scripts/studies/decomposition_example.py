"""One typical CaloParticle's decomposition into SimClusters, drawn once per collection.

    python scripts/studies/decomposition_example.py bound_electron_fc0_20GeV_testNanoML.root
    python scripts/studies/decomposition_example.py bound_photon_fc0_20GeV_testNanoML.root \
        --pattern g --fine-pattern "e+ e-"          # crosses Standard's surface, converts after it
    python scripts/studies/decomposition_example.py <fc0> <fc1> --event 124 --cp 1

What matters is whether a particle was created **before or after the boundary**:

- created after it (inside the volume, below a particle that crossed): it must
  not get a SimCluster of its own -- its energy is already in the crossing
  particle's. Giving it none is right.
- created before it, so it crossed the boundary itself: it should get its own
  SimCluster. Giving it none is wrong -- and its energy is then recorded nowhere.

The boundary collection opens a SimCluster only at a crossing track with no
crossing ancestor; everything below it gets no SimCluster of its own, its hits
go to the ancestor's and its energy is not recorded. It is therefore
right for particles created after the boundary and wrong for radiation emitted
before it, which crossed but has a crossing ancestor. The legacy collection
gives every flagged track its own SimCluster: right for every crossing particle,
but in Standard also wrong the other way -- tracks created inside CALO are
flagged at their vertex and get one (`interior`). FineCalo never flags a
track born inside.

Before or after is read from nanoML: a legacy SimCluster on a face of the
envelope crossed it (created before); one labelled `interior` was created
inside (after). Against Geant4's history from RECO (`check_boundary_nesting.py`,
γ 20 GeV) this split reproduces the energy without a boundary SimCluster per
event: created inside 20.62 GeV both ways (Standard), crossed 6.54 vs 6.71 GeV
(Standard) and 17.82 vs 18.07 GeV (FineCalo) -- the difference is crossings that
left no hits, which nanoML does not hold.

Takes one pair (`samples.py`) and needs its CaloParticle match
(`calopart_matching_<tag>.parquet`, the match step of scripts/analyze.py).

1. **Census** (printed): per config, the energy per event in roots, in crossings
   without a boundary SimCluster (wrong) and in particles created inside without
   one (right); and the boundary-collection pattern of each CaloParticle with
   its share.
2. **Pick the example**: by default the most common Standard pattern among
   CaloParticles with a crossing that has no boundary SimCluster, required in
   both configs (`--pattern`, `--fine-pattern` choose others); among those with
   at most the median number of legacy SimClusters, the one whose recorded share
   is nearest the median in both configs. `--event/--cp` (Standard row) override.
3. **Worked example** (printed): every SimCluster of that CaloParticle on both
   collections and configs, marked root / crossed / created inside, with the
   energy not recorded and the check that hits are kept.
4. **Figures**, each drawn twice on the same CaloParticle -- once per collection,
   same layout and axes, so the pair differs only in the collection:
   `decomposition_energy` (stacked bars of the energy crossing the boundary) and
   `decomposition_sketch` (schematic: before / after the boundary, and which
   particles get a SimCluster). Saved to outputs/<particle>/<energy>/<table>/
   (`_bound` for the boundary collection; a chosen pattern adds its name).

Energies are kinetic (`boundaryEnergy` minus rest mass): at the crossing, or for a particle
created inside at its creation. Which particle radiated which is not in the
file (no parent branches), so the sketch's branching points are drawn, not
measured; a particle without a SimCluster of its own is drawn from the root
whose impact point is nearest.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

from boundary_simclusters import samples
from boundary_simclusters.hep_style import (BASE_FONT, CONFIG_COLLECTION, GUIDE, INK, INK_2, MUTED,
                                            energy_text, save, series_label, single_figure,
                                            species, use_house_style)
from boundary_simclusters.output_paths import table_dir

KEYS, LABELS = samples.KEYS, samples.LABELS
COLLECTIONS = ("legacy", "boundary")
SUFFIX = {"legacy": "", "boundary": "_bound"}
PARTICLE_TEX = {"electron": "$e^-$", "photon": r"$\gamma$"}
SURFACE = {"fc0": "|z| = 293.5 cm", "fc1": "|z| = 321.05 cm"}
SHORT = {11: "e-", -11: "e+", 22: "g"}
SLUG = {"e-": "em", "e+": "ep", "g": "g"}
NONROOT_MIN = 0.01         # of the gun energy: a CaloParticle "has particles without a SimCluster" above this
SHOWN = 5                  # crossings without a SimCluster, drawn one by one; the softer ones are grouped
SHOWN_INSIDE = 3           # likewise for particles created inside
SUMMARY = "summary"        # gid of the text right of the sketch
LIGHT = {CONFIG_COLLECTION[("fc1", "boundary")]}      # gold: dark type on it


# --- reading ---------------------------------------------------------------------


def pattern(pdg):
    """The species of one CaloParticle's SimClusters: "e-", "e+ e-", "g x2" ..."""
    counts = pd.Series([SHORT.get(int(p), str(int(p))) for p in pdg]).value_counts()
    return " ".join(k if n == 1 else f"{k} x{n}" for k, n in sorted(counts.items()))


def with_roles(frames, key):
    """The legacy SimClusters with their role: `root` (also a boundary SimCluster),
    `crossed` (crossed the boundary, below a crossing ancestor) or `inside`
    (created inside the volume, below a crossing ancestor)."""
    bound, legacy = frames[key]["boundary"], frames[key]["legacy"]
    roots = set(zip(bound.event, bound.root_track))
    is_root = np.array([k in roots for k in zip(legacy.event, legacy.root_track)], dtype=bool)
    role = np.where(is_root, "root", np.where(legacy.born_inside, "inside", "crossed"))
    return legacy.assign(role=role)


def per_caloparticle(frames, legacy, key, e_gun):
    """One row per CaloParticle with SimClusters on both collections."""
    gb = frames[key]["boundary"].groupby(["event", "cp"])
    energy = legacy.pivot_table(index=["event", "cp"], columns="role", values="ekin",
                                aggfunc="sum", fill_value=0.0)
    energy = energy.reindex(columns=["root", "crossed", "inside"], fill_value=0.0)
    out = pd.DataFrame({"pattern": gb.pdg.agg(pattern), "recorded": gb.ekin.sum(),
                        "n_legacy": legacy.groupby(["event", "cp"]).size()}).dropna()
    out = out.join(energy[["crossed", "inside"]], how="left").fillna({"crossed": 0.0, "inside": 0.0})
    out["frac"] = out.recorded / e_gun
    out["lost"] = out.crossed / e_gun
    out["inside_frac"] = out.inside / e_gun
    return out


def gun_energy(sample):
    return float(sample.energy.removesuffix("GeV").replace("p", "."))


# --- census and pick ----------------------------------------------------------------


def census(cps, legacies, frames, sample):
    print(f"\n=== census: {sample.particle} {energy_text(sample.energy)} ({sample.tag})")
    print("energy per event [GeV]       roots (own SimCluster)   crossed, no SimCluster (wrong)   "
          "created inside, no SimCluster (right)")
    for key in KEYS:
        n = sample.tree(key).num_entries
        e = legacies[key].groupby("role").ekin.sum().reindex(["crossed", "inside"], fill_value=0) / n
        c = legacies[key].role.value_counts().reindex(["root", "crossed", "inside"], fill_value=0)
        roots = frames[key]["boundary"]
        print(f"  {LABELS[key]:9s} {roots.ekin.sum() / n:24.2f} {e['crossed']:20.2f} ({c['crossed']:6d} SC) "
              f"{e['inside']:20.2f} ({c['inside']:5d} SC)")
        print(f"  {'':9s} {len(roots) - c['root']} of {len(roots)} roots crossed but have no legacy "
              f"SimCluster")
    for key in KEYS:
        c = cps[key]
        table = (c.groupby("pattern")
                  .agg(n=("frac", "size"), recorded=("frac", "median"), lost=("lost", "median"),
                       inside=("inside_frac", "median"), n_legacy=("n_legacy", "median"))
                  .sort_values("n", ascending=False).head(5))
        table.insert(1, "share", table.n / len(c))
        print(f"\n{LABELS[key]}: {len(c)} CaloParticles, boundary-collection pattern "
              f"(medians as a share of the gun energy)")
        print(table.rename(columns={"recorded": "recorded", "lost": "crossed, no SC",
                                    "inside": "inside, no SC", "n_legacy": "legacy SCs"})
                   .to_string(float_format=lambda v: f"{v:.3f}"))


def matched(cps, sample):
    _, cm = sample.mappings()
    left = cps["fc0"].add_suffix("_0").rename_axis(["fc0_row", "fc0_cp"]).reset_index()
    right = cps["fc1"].add_suffix("_1").rename_axis(["fc1_row", "fc1_cp"]).reset_index()
    return cm.merge(left, on=["fc0_row", "fc0_cp"]).merge(right, on=["fc1_row", "fc1_cp"])


def percentile(population, value):
    return float((np.asarray(population) < value).mean())


def pick(cps, pairs, std_pattern=None, fine_pattern=None, event=None, cp=None):
    """The matched row of the example."""
    if event is not None:
        row = pairs[(pairs.fc0_row == event) & (pairs.fc0_cp == cp)]
        if row.empty:
            raise SystemExit(f"Standard event {event}, CaloParticle {cp} has no match in FineCalo")
        return row.iloc[0]
    std = cps["fc0"]
    std_pattern = std_pattern or std[std.lost > NONROOT_MIN].pattern.mode().iloc[0]
    fine_pattern = fine_pattern or std_pattern
    nonroot = lambda i: pairs[f"lost_{i}"] + pairs[f"inside_frac_{i}"] > NONROOT_MIN
    pool = pairs[(pairs.pattern_0 == std_pattern) & (pairs.pattern_1 == fine_pattern)
                 & nonroot(0) & nonroot(1)]
    if pool.empty:
        raise SystemExit(f"no CaloParticle with patterns {std_pattern!r} / {fine_pattern!r} and particles without a SimCluster")
    cand = pool[pool.n_legacy_0 <= pool.n_legacy_0.median()]
    score = sum((cand[f"frac_{i}"].map(lambda v, p=pool[f"frac_{i}"]: percentile(p, v)) - 0.5).abs()
                for i in (0, 1))
    row = cand.loc[score.idxmin()]
    print(f"\n=== pick: Standard pattern {std_pattern!r}, FineCalo pattern {fine_pattern!r}, with "
          f"particles without a SimCluster in both: {len(pool)} of {len(pairs)} matched CaloParticles ({len(pool) / len(pairs):.1%})")
    for i, key in ((0, "fc0"), (1, "fc1")):
        print(f"  {LABELS[key]}: recorded {row[f'frac_{i}']:.3f} of the gun energy, percentile "
              f"{percentile(pool[f'frac_{i}'], row[f'frac_{i}']):.2f} of the pool; "
              f"{int(row[f'n_legacy_{i}'])} legacy SimClusters (pool median {pool[f'n_legacy_{i}'].median():.0f})")
    return row


# --- the example ------------------------------------------------------------------


def example(frames, legacies, row):
    """Per config: the roots (boundary SimClusters, with whether the root also has
    a legacy SimCluster -- not, e.g., a photon that converts after crossing, whose
    hits go to the pair's legacy SimClusters), the crossings without a boundary SimCluster and
    the particles created inside, with their measured numbers."""
    out = {}
    for i, key in ((0, "fc0"), (1, "fc1")):
        ev, cp = int(row[f"fc{i}_row"]), int(row[f"fc{i}_cp"])
        sel = lambda f: f[(f.event == ev) & (f.cp == cp)]
        bound, legacy = sel(frames[key]["boundary"]), sel(legacies[key])
        roots = (bound.assign(in_legacy=bound.root_track.isin(legacy.root_track).to_numpy())
                      .sort_values("ekin", ascending=False))
        below = legacy[legacy.role != "root"].sort_values("ekin", ascending=False)
        # A particle without a boundary SimCluster is drawn from the root whose impact
        # point is nearest. Drawing only.
        d = np.hypot(below.x.to_numpy()[:, None] - roots.x.to_numpy()[None],
                     below.y.to_numpy()[:, None] - roots.y.to_numpy()[None])
        below = below.assign(owner=roots.root_track.to_numpy()[d.argmin(axis=1)] if len(roots) else -1)
        out[key] = dict(event=ev, cp=cp, bound=roots, legacy=legacy, roots=roots,
                        crossed=below[below.role == "crossed"], inside=below[below.role == "inside"])
    return out


def report(ex, sample, e_gun):
    print(f"\n=== worked example: {sample.particle} {energy_text(sample.energy)}, "
          f"Standard event {ex['fc0']['event']} CaloParticle {ex['fc0']['cp']}, "
          f"FineCalo event {ex['fc1']['event']} CaloParticle {ex['fc1']['cp']}")
    for key in KEYS:
        d = ex[key]
        print(f"\n{LABELS[key]}, legacy collection: {len(d['legacy'])} SimClusters "
              f"(root: also a boundary SimCluster; crossed: crossed the boundary, no boundary "
              f"SimCluster; inside: created inside the volume)")
        show = d["legacy"].sort_values(["role", "ekin"], ascending=[False, False])
        print(show[["role", "pdg", "ekin", "rec", "nhits", "face", "abs_z", "r"]]
              .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        b = d["bound"]
        print(f"{LABELS[key]}, boundary collection: {len(b)} SimCluster(s)")
        print(b[["pdg", "ekin", "rec", "nhits", "face", "in_legacy"]]
              .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        recorded, crossed, inside = b.ekin.sum(), d["crossed"].ekin.sum(), d["inside"].ekin.sum()
        print(f"  boundary collection: {len(b)} crossing particles with a SimCluster ({recorded:.3f} GeV); "
              f"{len(d['crossed'])} crossed without one (wrong) -> not recorded {crossed:.3f} GeV "
              f"({crossed / e_gun:.1%} of {e_gun:g} GeV); {len(d['inside'])} created inside, "
              f"no SimCluster (right, {inside:.3f} GeV, already in the crossing particle's)")
        missing = b[~b.in_legacy]
        print(f"  legacy collection: {len(d['legacy'])} SimClusters; {len(missing)} crossed without one "
              f"{' (wrong)' if len(missing) else ''}; {len(d['inside'])} for "
              f"particles created inside{' (wrong)' if len(d['inside']) else ''}")
        rec_b, rec_l = b.rec.sum(), d["legacy"].rec.sum()
        print(f"  hit energy: legacy {rec_l:.3f} GeV, boundary {rec_b:.3f} GeV -> hits kept: "
              f"{'yes' if np.isclose(rec_b, rec_l, rtol=1e-3) else 'NO'}")
        legacy_roots = d["legacy"][d["legacy"].role == "root"].set_index("root_track").ekin
        for r in b[b.in_legacy].itertuples():
            if not np.isclose(legacy_roots.loc[r.root_track], r.ekin, rtol=1e-4):
                print(f"  WARNING: root {r.root_track}: boundary {r.ekin:.3f} "
                      f"vs legacy {legacy_roots.loc[r.root_track]:.3f} GeV")


# --- figures ------------------------------------------------------------------------


def tex(pdg):
    return species(pdg)[0]


def identity(sample, ex):
    return (f"{PARTICLE_TEX[sample.particle]} {energy_text(sample.energy)}: Standard event "
            f"{ex['fc0']['event']}, FineCalo event {ex['fc1']['event']}, CaloParticle "
            f"{ex['fc0']['cp']}/{ex['fc1']['cp']}")


def segments(d, table):
    """(width, kind, label) along one bar -- the energy crossing the boundary:
    the roots, then the crossings without a boundary SimCluster (legacy: one segment each; boundary:
    one block not recorded). Particles created inside are not stacked: their energy
    is part of the crossing particle's."""
    segs = [(r.ekin, "root" if table == "boundary" or r.in_legacy else "lost",
             f"{tex(r.pdg)} {r.ekin:.1f}" if table == "boundary" or r.in_legacy
             else f"{tex(r.pdg)} {r.ekin:.1f}: no SimCluster")
            for r in d["roots"].itertuples()]
    crossed = d["crossed"]
    if table == "boundary":
        if len(crossed):
            segs.append((crossed.ekin.sum(), "lost", "crossed, no SimCluster"))
        return segs
    for r in crossed.head(SHOWN).itertuples():
        segs.append((r.ekin, "crossed", f"{tex(r.pdg)} {r.ekin:.1f}"))
    rest = crossed.iloc[SHOWN:]
    if len(rest):
        segs.append((rest.ekin.sum(), "crossed", f"+{len(rest)}"))
    return segs


def figure_energy(ex, sample, table, e_gun, title_id, xmax):
    fig, ax = single_figure(figsize=(10.0, 5.8))
    height, handles = 0.5, []
    for y, key in ((1, "fc0"), (0, "fc1")):
        d, colour = ex[key], CONFIG_COLLECTION[(key, table)]
        text_colour = INK if colour in LIGHT else "white"
        left = 0.0
        for width, kind, label in segments(d, table):
            if kind == "lost":
                ax.barh(y, width, left=left, height=height, facecolor="none", edgecolor=colour,
                        hatch="///", linewidth=1.4)
                if width > 0.2 * xmax:
                    ax.text(left + width / 2, y, label, ha="center", va="center",
                            fontsize=BASE_FONT - 3, color=INK,
                            bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
            else:
                ax.barh(y, width, left=left, height=height, color=colour, edgecolor="white",
                        linewidth=1.5)
                if width > 0.07 * xmax:
                    ax.text(left + width / 2, y, label, ha="center", va="center",
                            fontsize=BASE_FONT - 3, color=text_colour)
            left += width
        ax.plot(d["bound"].rec.sum(), y + height / 2 + 0.11, marker="v", color=INK, markersize=9,
                linestyle="none")
        crossing = d["roots"].ekin.sum() + d["crossed"].ekin.sum()
        if table == "legacy":
            numbers = (f"{len(d['legacy'])} SimClusters; crossing the boundary: "
                       f"Σ$E_\\mathrm{{kin}}$ {crossing:.1f} GeV")
        else:
            numbers = (f"{len(d['bound'])} SimCluster{'s' * (len(d['bound']) > 1)}, "
                       f"$E_\\mathrm{{kin}}$ {d['bound'].ekin.sum():.1f} GeV recorded")
        handles.append(Patch(facecolor=colour, edgecolor="white",
                             label=f"{series_label(key, table, sample.energy)}  ({numbers})"))
    handles.append(Line2D([], [], color=INK, marker="v", markersize=9, linestyle="none",
                          label=f"hit energy (Σ recEnergy): Standard {ex['fc0']['bound'].rec.sum():.1f}, "
                                f"FineCalo {ex['fc1']['bound'].rec.sum():.1f} GeV"))
    inside = {k: ex[k]["inside"] for k in KEYS if len(ex[k]["inside"])}
    if table == "boundary":
        lost = {k: ex[k]["crossed"].ekin.sum() for k in KEYS}
        handles.append(Patch(facecolor="none", edgecolor=INK_2, hatch="///",
                             label=f"crossed the boundary, no SimCluster of its own: not recorded "
                                   f"(Standard {lost['fc0']:.1f}, FineCalo {lost['fc1']:.1f} GeV)"))
        for k, f in inside.items():
            handles.append(Patch(facecolor="none", edgecolor="none",
                                 label=f"{LABELS[k]}: {len(f)} created after the boundary, without a SimCluster "
                                       f"({f.ekin.sum():.1f} GeV, already in the crossing particle's)"))
    else:
        missing = {k: ex[k]["roots"][~ex[k]["roots"].in_legacy] for k in KEYS}
        if any(len(m) for m in missing.values()):
            handles.append(Patch(facecolor="none", edgecolor=INK_2, hatch="///",
                                 label=f"crossed the boundary, no SimCluster of its own "
                                       f"(its hits are in other SimClusters)"))
        else:
            handles.append(Patch(facecolor="none", edgecolor="none",
                                 label="every particle crossing the boundary has its own SimCluster"))
        for k, f in inside.items():
            handles.append(Patch(facecolor="none", edgecolor="none",
                                 label=f"{LABELS[k]}: {len(f)} more for particles created after "
                                       f"the boundary ({f.ekin.sum():.1f} GeV, not stacked)"))
    ax.axvline(e_gun, ymax=0.55, color=INK_2, linewidth=1.3, linestyle=":")   # stops below the legend
    ax.annotate("gun energy", xy=(e_gun, 0.02), xycoords=("data", "axes fraction"), xytext=(-5, 0),
                textcoords="offset points", ha="right", va="bottom", fontsize=BASE_FONT - 4, color=INK_2)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.5, 2.75)               # fixed, so the two collections' figures line up
    ax.set_yticks([1, 0], ["Standard", "FineCalo"])
    ax.tick_params(axis="y", which="both", length=0)
    ax.minorticks_on()
    ax.tick_params(axis="y", which="minor", left=False, right=False)
    ax.legend(handles=handles, loc="upper left", fontsize=BASE_FONT - 4)
    ax.set_xlabel("kinetic energy at the crossing  [GeV]\n(one segment per SimCluster)",
                  linespacing=1.5)
    ax.set_title(f"{table.capitalize()} collection: one CaloParticle\n{title_id}",
                 pad=12, fontsize=BASE_FONT)
    return fig


def lane(ax, d, key, table, yc, half, e_gun):
    """One config's view of the CaloParticle, schematically: particles created
    before the boundary (left of it) and after it (inside, the shaded band)."""
    colour = CONFIG_COLLECTION[(key, table)]
    x_surf, x_conv, x_in, x_label, x_box, x_group = 6.0, 6.45, 7.4, 7.65, 12.7, 13.1
    roots, crossed, inside = d["roots"], d["crossed"], d["inside"]
    width = lambda e: 1.0 + 3.0 * np.sqrt(max(e, 0.0) / e_gun)
    style = lambda pdg: ":" if pdg == 22 else "-"
    c_shown, c_rest = crossed.head(SHOWN), crossed.iloc[SHOWN:]
    i_shown, i_rest = inside.head(SHOWN_INSIDE), inside.iloc[SHOWN_INSIDE:]

    # Endpoints top to bottom: each root, the crossings drawn from it, then what
    # was created after it crossed. A gap separates one root's group from the next.
    order, gaps = [], []
    for r in roots.itertuples():
        mine = lambda f: f[f.owner == r.root_track]
        order.append(("root", r))
        order += [("crossed", f) for f in mine(c_shown).itertuples()]
        if len(mine(c_rest)):
            order.append(("crossed_rest", mine(c_rest)))
        order += [("inside", f) for f in mine(i_shown).itertuples()]
        if len(mine(i_rest)):
            order.append(("inside_rest", mine(i_rest)))
        gaps.append(len(order))
    slots = np.arange(len(order)) + 1.2 * np.searchsorted(gaps, np.arange(len(order)), side="right")
    span = slots[-1] if len(order) > 1 else 1.0
    ys = yc + half - 2 * half * slots / span if len(order) > 1 else np.array([yc])
    end = dict(enumerate(ys))
    step = 2 * half / span
    owner_of = lambda kind, f: f.owner.iloc[0] if kind.endswith("rest") else f.owner
    energy_of = lambda kind, f: f.ekin.sum() if kind.endswith("rest") else f.ekin
    pdg_of = lambda kind, f: 22 if kind.endswith("rest") else f.pdg

    # Before / after the boundary.
    top, bottom = yc + half + 0.22, yc - half - 0.22
    ax.fill_between([x_surf, x_in + 0.05], bottom, top, color=GUIDE, alpha=0.3, linewidth=0)
    ax.plot([x_surf, x_surf], [bottom - 0.05, top + 0.05], color=INK, linewidth=1.6)
    ax.text(x_surf, top + 0.08, f"{LABELS[key]}'s surface, {SURFACE[key]}", ha="center",
            va="bottom", fontsize=BASE_FONT - 4)
    ax.text(x_surf - 0.15, bottom - 0.03, "created before the boundary", ha="right", va="top",
            fontsize=BASE_FONT - 5, color=INK_2)
    ax.text((x_surf + x_in) / 2, bottom - 0.03, "after", ha="center", va="top",
            fontsize=BASE_FONT - 5, color=INK_2)
    ax.plot(0, yc, marker="o", color=INK, markersize=5)
    start = (0.0, yc)
    if len(roots) > 1:                     # the gun photon converts before the surface
        start = (1.4, yc)
        ax.plot([0, 1.4], [yc, yc], color=colour, linestyle=":", linewidth=width(e_gun))
        ax.text(0.7, yc + 0.06, r"$\gamma$ (no crossing)", ha="center", va="bottom", fontsize=BASE_FONT - 5)

    # Root lines; crossings branch off before the surface (earliest -> lowest
    # endpoint), particles created inside off the root's extension past it.
    for i, (kind, r) in enumerate(order):
        if kind != "root":
            continue
        p0, p1 = start, (x_surf, end[i])
        ax.plot([p0[0], x_surf], [p0[1], end[i]], color=colour, linestyle=style(r.pdg),
                linewidth=width(r.ekin), solid_capstyle="butt")
        before = [j for j, (k, f) in enumerate(order)
                  if k.startswith("crossed") and owner_of(k, f) == r.root_track]
        after = [j for j, (k, f) in enumerate(order)
                 if k.startswith("inside") and owner_of(k, f) == r.root_track]
        xs = np.linspace(p0[0] + 0.8, x_surf - 1.2, len(before) + 2)[1:-1]
        for x_b, j in zip(xs, reversed(before)):
            y_b = p0[1] + (p1[1] - p0[1]) * (x_b - p0[0]) / (p1[0] - p0[0])
            kind, f = order[j]
            ax.plot([x_b, x_surf], [y_b, end[j]], color=colour, linestyle=style(pdg_of(kind, f)),
                    linewidth=width(energy_of(kind, f)) if kind == "crossed" else 1.0)
        if after:
            ax.plot([x_surf, x_conv], [end[i], end[i]], color=colour, linestyle=style(r.pdg),
                    linewidth=width(r.ekin))
            for j in after:
                kind, f = order[j]
                ax.plot([x_conv, x_in], [end[i], end[j]], color=colour, linestyle=style(pdg_of(kind, f)),
                        linewidth=width(energy_of(kind, f)) if kind == "inside" else 1.0)

    # Labels: identical in both collections, except that a particle this
    # collection gives no SimCluster of its own is set light and italic.
    for i, (kind, f) in enumerate(order):
        if kind.endswith("rest"):
            text = f"+{len(f)} softer, Σ {f.ekin.sum():.2f} GeV"
        else:
            text = f"{tex(f.pdg)}  {f.ekin:.2f} GeV"
        if kind.startswith("inside"):
            text += ", created inside"
        no_legacy = kind == "root" and not f.in_legacy
        if table == "legacy" and no_legacy:
            text += ": no SimCluster"
        own = kind == "root" and f.in_legacy if table == "legacy" else kind == "root"
        own = own or table == "legacy" and kind != "root"
        ax.text(x_label, end[i], text, va="center", fontsize=BASE_FONT - 4,
                color=INK if own else MUTED, fontstyle="normal" if own else "italic",
                fontweight="bold" if kind == "root" and own else None)
        if table == "legacy" and not no_legacy:
            x_mark = x_in if kind.startswith("inside") else x_surf
            ax.plot(x_mark, end[i], marker="D" if kind.startswith("inside") else "s", color=colour,
                    markersize=6.5, markeredgecolor=INK, markeredgewidth=0.6)

    # What this collection made of the CaloParticle, in the same place in both versions.
    created = len(d["legacy"]) if table == "legacy" else len(roots)
    missing = int((~roots.in_legacy).sum()) if table == "legacy" else len(crossed)
    lines = [(f"{created} SimCluster{'s' * (created != 1)} created", True),
             (f"{missing} crossed the boundary without a SimCluster", bool(missing))]
    if table == "legacy":
        lines.append((f"{len(inside)} SimCluster{'s' * (len(inside) != 1)} "
                       f"for particles created after it", bool(len(inside))))
        lines.append((f"$E_\\mathrm{{kin}}$ recorded: {d['legacy'].ekin.sum():.2f} GeV", False))
    else:
        lines.append((f"{len(inside)} created after it without a SimCluster", False))
        lines.append((f"$E_\\mathrm{{kin}}$ recorded: {roots.ekin.sum():.2f} GeV, "
                      f"not recorded: {crossed.ekin.sum():.2f} GeV", False))
    for k, (text, bold) in enumerate(lines):
        ax.text(x_group, yc + 0.33 - 0.22 * k, text, va="center", fontsize=BASE_FONT - 3,
                fontweight="bold" if bold or k == 0 else None, gid=SUMMARY)
    if table == "legacy":
        return

    # Boundary: one box per root, around everything below it that has no SimCluster of its own.
    groups = [i for i, (k, _) in enumerate(order) if k == "root"] + [len(order)]
    pad = min(0.52 * step, 0.15)
    for a, b in zip(groups[:-1], groups[1:]):
        top_b, bottom_b = end[a] + pad, end[b - 1] - pad
        ax.add_patch(FancyBboxPatch((x_surf - 0.15, bottom_b), x_box - (x_surf - 0.15),
                                    top_b - bottom_b, boxstyle="round,pad=0.0,rounding_size=0.06",
                                    facecolor="none", edgecolor=colour, linewidth=2.0))


def figure_sketch(ex, sample, table, e_gun, title_id):
    fig, ax = single_figure(figsize=(14.2, 7.8))
    ax.axis("off")
    half = 0.72
    for yc, key in ((2.75, "fc0"), (0.2, "fc1")):
        lane(ax, ex[key], key, table, yc, half, e_gun)
    ax.set_xlim(-0.3, 18.0)
    ax.set_ylim(-1.05, 4.0)
    ax.set_title(f"{table.capitalize()} collection: one CaloParticle\n{title_id}", pad=8,
                 fontsize=BASE_FONT)
    return fig


def draw(ex, sample, e_gun, slug):
    title_id = identity(sample, ex)
    xmax = 1.08 * max([e_gun] + [ex[k]["roots"].ekin.sum() + ex[k]["crossed"].ekin.sum() for k in KEYS]
                      + [ex[k]["bound"].rec.sum() for k in KEYS])
    tree = sample.tree("fc0")
    written = []
    for table in COLLECTIONS:
        out = table_dir(tree, table)
        name = f"{slug}{sample.tag}{SUFFIX[table]}.png"
        written.append(save(figure_energy(ex, sample, table, e_gun, title_id, xmax),
                            out / f"decomposition_energy_{name}"))
        written.append(save(figure_sketch(ex, sample, table, e_gun, title_id),
                            out / f"decomposition_sketch_{name}"))
    return written


def pattern_slug(std, fine):
    """"g", "e+ e-" -> "g-epem_": names the figures of a chosen pattern."""
    short = lambda p: "".join(SLUG.get(s, s) for s in p.split() if not s.startswith("x"))
    return f"{short(std)}-{short(fine or std)}_"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", type=Path, help="one pair of nanoML files (or one; its partner is found)")
    ap.add_argument("--tag", help="sample tag (default: the match step's)")
    ap.add_argument("--pattern", help='Standard boundary-collection pattern to pick from, e.g. "g", "e+ e-"')
    ap.add_argument("--fine-pattern", help="FineCalo's pattern for the same CaloParticle (default: --pattern)")
    ap.add_argument("--event", type=int, help="Standard row of the CaloParticle to draw (default: picked)")
    ap.add_argument("--cp", type=int, default=0, help="its CaloParticle index in Standard (with --event)")
    args = ap.parse_args()
    if len(args.files) > 2:
        raise SystemExit("one pair at a time")
    sample = samples.load_sample(*args.files, tag=args.tag)
    if "boundary" not in sample.tables:
        raise SystemExit(f"{sample.tag}: no boundary collection in these files")
    e_gun = gun_energy(sample)
    frames = {key: {t: samples.simcluster_frame(sample.tree(key), t, key) for t in COLLECTIONS}
              for key in KEYS}
    legacies = {key: with_roles(frames, key) for key in KEYS}
    cps = {key: per_caloparticle(frames, legacies[key], key, e_gun) for key in KEYS}
    census(cps, legacies, frames, sample)
    row = pick(cps, matched(cps, sample), args.pattern, args.fine_pattern, args.event, args.cp)
    ex = example(frames, legacies, row)
    report(ex, sample, e_gun)
    use_house_style()
    slug = pattern_slug(args.pattern, args.fine_pattern) if args.pattern else ""
    print()
    for path in draw(ex, sample, e_gun, slug):
        print("wrote", path)


if __name__ == "__main__":
    main()
