"""The standard run's figures for one fc0/fc1 pair: file in, figures out.

    python scripts/analyze.py bound_photon_fc0_200GeV_testNanoML.root --steps figures
    python scripts/analyze.py a.root b.root --steps figures --table legacy --tag s1
    python scripts/analyze.py a.root --steps figures --only energy n_per_event

`run(sample)` is the figures step of scripts/analyze.py.

Give one file and its partner is found (`samples.py`); give two in either order.
Config, particle, energy and tables are read from the files, never their names.
The boundary table (`BoundarySimCluster_*`) is read when both files carry it;
`--table legacy` reads `SimCluster_*`. Axis
labels name the table's own branches.

Every figure lands in one folder per sample and table,
outputs/<particle>/<energy>/<table>/<name>_<tag>.png (`output_paths.table_dir`):

    position_z        |z| of every crossing: each config's own front face
    position_z_off    the crossings off that face (bore and wall)
    eta               impactPoint_eta, the momentum direction (docs/data_handling.md)
    energy, pt        boundaryEnergy and boundaryPt, log-log, medians marked
    type_count        particle-type mix of the SimClusters
    type_momentum     particle-type share of the momentum crossing
    type_ratio        FineCalo/Standard count above a momentum threshold, per type
    n_per_event       SimClusters per event
    n_per_caloparticle SimClusters per CaloParticle (every CaloParticle of the file)
    ekin_per_caloparticle  summed SimCluster kinetic energy per CaloParticle over the gun
                      energy: ~1 on the legacy table, well below on the boundary table
                      (the energy of folded-in particles is not booked)
    rec_over_ekin     hit energy / kinetic energy at the crossing, per SimCluster
                      above 100 MeV: ~0.95 on the legacy table, far above 1 on the
                      boundary table, whose SimClusters hold hits of particles whose
                      energy they do not carry
    rec_energy        summed SimCluster recEnergy per CaloParticle (exclusive deposits,
                      so the sum is valid -- docs/data_handling.md)
    energy_weighted_fraction, energy_weighted_gev
                      kinetic energy crossing, weighted by E_kin

One look for every figure (`hep_style`): a single panel, step lines, colour by
config and table (`config_colours`: purple Standard / gold FineCalo on the
boundary table, blue / red on the legacy one), no experiment label, and per-series numbers in the
legend entry rather than drawn into the plot. The tag in the
filename is `--tag` or the sample's (`samples.py`), plus `_bound` on the
boundary table.

`boundaryPmag` is not drawn as a distribution: E = |p| for 99.84-99.90% of
SimClusters, so it would repeat `energy`. The energy-weighted spectrum uses
kinetic energy, not `boundaryEnergy`, which includes rest mass
(`boundary_geometry.kinetic_energy`).
"""


import numpy as np

from boundary_simclusters import samples  # noqa: E402
from boundary_simclusters.boundary_geometry import front_face, front_face_mask  # noqa: E402
from boundary_simclusters.hep_style import (BASE_FONT, CONFIG_COLLECTION, GUIDE, INK_2, MUTED,  # noqa: E402
                                            config_colours, energy_linestyle, series_label,
                       legend_with_headroom, save, single_figure, use_house_style)
from boundary_simclusters.output_paths import style_dir, table_dir  # noqa: E402

# FineCalo first, as in every legend of the repo; colour is bound to config and table.
# A series key is a config ("fc1") in a figure of one table, or (config, table)
# in a comparison figure, which overlays both tables (`compare`).
KEYS = ("fc1", "fc0")
COMPARE_KEYS = (("fc1", "boundary"), ("fc0", "boundary"), ("fc1", "legacy"), ("fc0", "legacy"))
LABELS = samples.LABELS
COLORS = config_colours("boundary")       # reset per run to the table drawn
SERIES = {}                               # legend names, set per run
LINESTYLE = "-"                           # the sample's energy, per run
TABLE_NOW, ENERGY_NOW = "boundary", ""    # for legends without a per-config series
# The font lacks U+207A/U+207B, so particle labels are mathtext.
PARTICLE_TEX = {"electron": "$e^-$", "positron": "$e^+$", "photon": r"$\gamma$"}
FONT = BASE_FONT

# Set once per run in `draw`: the folder figures go to, and the branch prefix
# of the table drawn, which every axis label names.
FIGDIR = None
P = "SimCluster"


def finish(fig, ax, name, *, title, xlabel, ylabel, loc="upper right"):
    """The parts every figure shares: labels, title, legend."""
    ax.set_xlabel(xlabel, linespacing=1.5)
    ax.set_ylabel(ylabel, linespacing=1.5)
    ax.set_title(title, pad=12)
    legend_with_headroom(ax, loc=loc)
    return save(fig, FIGDIR / name)


def steps(ax, values, bins, key, label, **kw):
    """One config's step histogram. Returns the counts."""
    counts, edges = np.histogram(values, bins=bins)
    kw.setdefault("linestyle", LINESTYLE)
    ax.stairs(counts, edges, color=COLORS[key], label=label, **kw)
    return counts


def config_of(key):
    """The config of a series key: "fc1", or the first part of ("fc1", "legacy")."""
    return key if isinstance(key, str) else key[0]


def who(key):
    """A series in running text: "Standard", or "Standard (legacy)"."""
    return LABELS[config_of(key)] + ("" if isinstance(key, str) else f" ({key[1]})")


def z_label():
    return "$|z|$ of impact point  [cm]"


# --- where the crossings sit --------------------------------------------------


def on_face(d):
    """(front face |z|, mask of crossings on it) for one config."""
    plane = front_face(d.abs_z.to_numpy())
    return plane, front_face_mask(d.abs_z.to_numpy(), plane)


def on_plane_title(frames):
    """Does each config actually record on a plane? Computed, never hard-coded:
    the photon pair puts only 86% of Standard's crossings on one."""
    shares = {k: float(on_face(frames[k])[1].mean()) for k in KEYS}
    worst = min(shares, key=shares.get)
    if shares[worst] > 0.95:
        return "Each config records on its own front face"
    return f"Only {100 * shares[worst]:.0f}% of {who(worst)}'s crossings are on its front face"


def figure_position_z(frames, tag):
    fig, ax = single_figure()
    pooled = np.concatenate([frames[k].abs_z for k in KEYS])
    bins = np.arange(np.floor(pooled.min()) - 1, np.ceil(pooled.max()) + 1, 0.5)
    for key in KEYS:
        plane, mask = on_face(frames[key])
        steps(ax, frames[key].abs_z, bins, key,
              f"{SERIES[key]}  ({100 * mask.mean():.2f}% at {plane:.1f} cm)")
    ax.set_yscale("log")
    return finish(fig, ax, f"position_z_{tag}.png", title=on_plane_title(frames),
                  xlabel=z_label(), ylabel="SimClusters")


def figure_position_z_off(frames, tag):
    """The crossings off each config's own front face. They are not off the
    boundary: they sit on the envelope's bore and wall faces."""
    fig, ax = single_figure()
    pooled = np.concatenate([frames[k].abs_z for k in KEYS])
    # 2 cm bins: at 1 cm the tail is mostly bins of one and reads as noise.
    bins = np.arange(np.floor(pooled.min()) - 1, np.ceil(pooled.max()) + 2, 2.0)
    for key in KEYS:
        plane, mask = on_face(frames[key])
        off = frames[key].abs_z[~mask]
        steps(ax, off, bins, key, f"{SERIES[key]}  ({len(off)} crossings)")
        # The face each tail departs from.
        ax.axvline(plane, color=GUIDE, linestyle=":", linewidth=1.2)
    return finish(fig, ax, f"position_z_off_{tag}.png",
                  title="Crossings off each config's front face",
                  xlabel=z_label(), ylabel="SimClusters")


# --- direction, energy, momentum -----------------------------------------------


def figure_eta(frames, tag):
    fig, ax = single_figure()
    bins = np.linspace(-6, 6, 121)
    for key in KEYS:
        abs_eta = frames[key].eta.abs()
        # max |eta| swings on one outlier, so the 99th percentile comes first.
        steps(ax, frames[key].eta, bins, key,
              f"{SERIES[key]}  (99%: |η| < {np.percentile(abs_eta, 99):.2f})")
    ax.set_yscale("log")
    return finish(fig, ax, f"eta_{tag}.png", title="Direction of the boundary crossings",
                  xlabel="$\\eta$ of the momentum at the crossing",
                  ylabel="SimClusters")


def log_steps(frames, field, lo, hi, fig_ax, nbins=60):
    """Log-x, log-y step histograms; each config's median is in its legend entry."""
    fig, ax = fig_ax
    bins = np.logspace(np.log10(lo), np.log10(hi), nbins)
    for key in KEYS:
        values = frames[key][field].to_numpy()
        values = values[values > 0]
        median = float(np.median(values))
        steps(ax, values, bins, key, f"{SERIES[key]}  (median {median:.3g} GeV)")
    ax.set_xscale("log")
    ax.set_yscale("log")


def figure_energy(frames, tag):
    fig, ax = single_figure()
    log_steps(frames, "energy", 1e-4, 300, (fig, ax))
    return finish(fig, ax, f"energy_{tag}.png", title="Energy carried across the boundary",
                  xlabel="$E$ at the crossing  [GeV]",
                  ylabel="SimClusters")


def figure_pt(frames, tag):
    fig, ax = single_figure()
    log_steps(frames, "pt", 1e-5, 100, (fig, ax))
    return finish(fig, ax, f"pt_{tag}.png", title="Transverse momentum at the crossing",
                  xlabel="$p_\\mathrm{T}$ at the crossing  [GeV]",
                  ylabel="SimClusters")


# --- particle type -----------------------------------------------------------------
# photon + e± is 99.8% of boundary SimClusters, so the share figures cover those
# three and give the hadronic remainder in the legend -- on a linear axis a
# 0.18% bar is invisible.

PARTICLE_TYPES = (("photon", 22), ("$e^-$", 11), ("$e^+$", -11))
# The ratio figure has no per-config series, so colour encodes particle type
# there -- deliberately outside the blue/orange config pair.
TYPE_COLOUR = {"photon": "#c99700", "$e^-$": "#199c8e", "$e^+$": "#7d4bb0"}
# Past this point Standard keeps too few objects for the ratio to mean anything.
MIN_COUNT = 25
# A mix this close counts as "the same": the electron sample's largest gap is 3.9 pp.
SAME_MIX_PP = 5.0


def shares(d):
    """Per-type share of objects and of summed momentum, one config."""
    pdg, p = d.pdg.to_numpy(), d.pmag.to_numpy()
    counts = {name: int((pdg == code).sum()) for name, code in PARTICLE_TYPES}
    return {
        "count_share": {n: 100 * c / len(pdg) for n, c in counts.items()},
        "momentum_share": {n: 100 * float(p[pdg == code].sum()) / p.sum()
                           for n, code in PARTICLE_TYPES},
        "other_share": 100 * (len(pdg) - sum(counts.values())) / len(pdg),
    }


def mix_title(s):
    if "fc1" not in s:                      # a comparison figure: no single pair to rank
        return "Particle-type mix at the crossing, both collections"
    diffs = {n: s["fc1"]["count_share"][n] - s["fc0"]["count_share"][n] for n, _ in PARTICLE_TYPES}
    biggest = max(diffs, key=lambda n: abs(diffs[n]))
    if abs(diffs[biggest]) < SAME_MIX_PP:
        return "Both configs cross the boundary with the same particle-type mix"
    direction = "more" if diffs[biggest] > 0 else "fewer"
    return f"The particle-type mix differs: {abs(diffs[biggest]):.0f} pp {direction} {biggest} in FineCalo"


def momentum_title(s):
    if "fc1" not in s:
        return "Momentum share per particle type, both collections"
    diffs = {n: s["fc1"]["momentum_share"][n] - s["fc0"]["momentum_share"][n] for n, _ in PARTICLE_TYPES}
    winner = max(diffs, key=diffs.get)
    if diffs[winner] < SAME_MIX_PP:
        return "Both configs share the momentum across particle types alike"
    return f"FineCalo moves {diffs[winner]:+.0f} pp of the momentum onto {winner}"


def type_bars(ax, s, field, labels):
    x = np.arange(len(PARTICLE_TYPES))
    width = 0.76 / len(KEYS)
    tallest = 0.0
    for i, key in enumerate(KEYS):
        values = [s[key][field][name] for name, _ in PARTICLE_TYPES]
        tallest = max(tallest, max(values))
        offset = (i - (len(KEYS) - 1) / 2) * width
        ax.bar(x + offset, values, width * 0.94, color=COLORS[key], label=labels[key])
        for xi, value in enumerate(values):
            ax.annotate(f"{value:.1f}%", xy=(xi + offset, value), xytext=(0, 4),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=FONT - 4 if len(KEYS) <= 2 else FONT - 6, color=INK_2)
    ax.set_xticks(x)
    ax.set_xticklabels([name for name, _ in PARTICLE_TYPES], fontsize=FONT - 1)
    ax.set_ylim(0, tallest * 1.10)


def figure_type_count(frames, tag):
    fig, ax = single_figure()
    s = {k: shares(frames[k]) for k in KEYS}
    type_bars(ax, s, "count_share",
              {k: f"{SERIES[k]}  (hadrons + other {s[k]['other_share']:.2f}%)" for k in KEYS})
    return finish(fig, ax, f"type_count_{tag}.png", title=mix_title(s),
                  xlabel="Particle type at the crossing",
                  ylabel="Share of SimClusters  [%]")


def figure_type_momentum(frames, tag):
    fig, ax = single_figure()
    s = {k: shares(frames[k]) for k in KEYS}
    type_bars(ax, s, "momentum_share", SERIES)
    return finish(fig, ax, f"type_momentum_{tag}.png", title=momentum_title(s),
                  xlabel="Particle type at the crossing",
                  ylabel="Share of $\\Sigma\\,|\\vec{p}\\,|$ at the crossing  [%]")


def ratio_pairs():
    """(legend suffix, FineCalo key, Standard key, line options) per ratio drawn:
    one in a figure of one table, one per table in a comparison figure, the legacy
    one marked with dots (line style is reserved for the energy)."""
    if isinstance(KEYS[0], str):
        return [("", "fc1", "fc0", {})]
    return [(f", {t}", ("fc1", t), ("fc0", t),
             {} if t == "boundary" else {"marker": "o", "markersize": 4, "markevery": 4})
            for t in ("boundary", "legacy")]


def figure_type_ratio(frames, tag):
    """How the FineCalo/Standard multiplicity ratio depends on momentum, per type."""
    fig, ax = single_figure()
    thresholds = np.logspace(-3, 2.1, 60)
    for suffix, fine, std, style in ratio_pairs():
        for name, code in PARTICLE_TYPES:
            p = {k: frames[k].pmag.to_numpy()[frames[k].pdg.to_numpy() == code] for k in (fine, std)}
            above = {k: (p[k][:, None] > thresholds[None, :]).sum(axis=0) for k in (fine, std)}
            # Past the point where Standard keeps only a few dozen objects the ratio
            # is counting noise; it is dropped rather than drawn.
            ratio = np.where(above[std] >= MIN_COUNT, above[fine] / np.maximum(above[std], 1), np.nan)
            ax.plot(thresholds, ratio, color=TYPE_COLOUR[name], linewidth=2.4,
                    linestyle=LINESTYLE, label=f"{name}{suffix or ', ' + TABLE_NOW}, {ENERGY_NOW}", **style)
        no_cut = len(frames[fine]) / len(frames[std])
        ax.axhline(no_cut, color=MUTED, linewidth=1.3, linestyle=":")
        # Anchored to the right edge of the axes: the curves end where Standard
        # runs out of objects, well before the last threshold.
        ax.annotate(f"all particle types, no cut{suffix}: {no_cut:.2f}×",
                    xy=(0.99, no_cut), xycoords=("axes fraction", "data"), xytext=(0, 5),
                    textcoords="offset points", ha="right", fontsize=FONT - 4, color=MUTED)
    ax.axhline(1.0, color=INK_2, linewidth=1.3, linestyle=(0, (1, 3)))
    ax.annotate("parity", xy=(0.99, 1.0), xycoords=("axes fraction", "data"), xytext=(0, -14),
                textcoords="offset points", ha="right", fontsize=FONT - 4, color=INK_2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # Six legend entries in a comparison figure: at the top, the axis grown to clear them.
    return finish(fig, ax, f"type_ratio_{tag}.png", loc="lower left" if len(KEYS) <= 2 else "upper left",
                  title="Object-count ratio above a momentum threshold, by particle type",
                  xlabel="Threshold on $|\\vec{p}\\,|$ at the crossing  [GeV]",
                  ylabel="FineCalo / Standard object count")


# --- multiplicity ----------------------------------------------------------------


def per_event(frame, n_events):
    """SimClusters per event, events with none included."""
    return np.bincount(frame.event.to_numpy(), minlength=n_events)


def count_bins(n_max):
    """Integer-centred bins, one wide unless that would exceed ~60 bins: the
    boundary table has 1-4 SimClusters per event, the legacy table ~50."""
    width = max(1, int(np.ceil((n_max + 1) / 60)))
    return np.arange(-0.5, n_max + width + 0.5, width)


def per_caloparticle(frame, tree, field=None):
    """Per CaloParticle, every CaloParticle of the file (those with no SimCluster
    included): the number of SimClusters, or the sum of `field` over them. Per
    file, like `per_event`: no Standard/FineCalo pairing."""
    n_cp = np.array([len(x) for x in tree["CaloPart_pdgId"].array(library="np")])
    first = np.concatenate([[0], np.cumsum(n_cp)[:-1]])
    event, cp = frame.event.to_numpy(), frame.cp.to_numpy()
    valid = (cp >= 0) & (cp < n_cp[event])
    weights = None if field is None else frame[field].to_numpy()[valid]
    return np.bincount(first[event[valid]] + cp[valid], weights=weights, minlength=int(n_cp.sum()))


def figure_ekin_per_caloparticle(ekin, e_gun, tag):
    """Summed kinetic energy of each CaloParticle's SimClusters over the gun
    energy: the energy each table books at the boundary. The legacy table closes
    (~1); the boundary table books only each first crossing track's energy, not
    that of the particles folded into it. Kinetic, since rest mass is not
    deposited; against the gun energy the choice moves the mean by < 0.001."""
    fig, ax = single_figure()
    top = max(1.1, 1.02 * max((ekin[k] / e_gun).max() for k in KEYS))
    # Bins centred on multiples of 0.02, so a CaloParticle booking the full gun
    # energy (an unconverted photon) falls in one bin centred on 1.
    bins = np.arange(-0.01, top + 0.02, 0.02)
    for key in KEYS:
        f = ekin[key] / e_gun
        steps(ax, f, bins, key, f"{SERIES[key]}  (mean {f.mean():.2f})")
    ax.set_xlim(0.0, bins[-1])
    ax.set_yscale("log")
    return finish(fig, ax, f"ekin_per_caloparticle_{tag}.png", loc="upper left",
                  title="Energy booked at the boundary per CaloParticle",
                  xlabel="$\\Sigma E_\\mathrm{kin}$ / $E_\\mathrm{gun}$ per CaloParticle",
                  ylabel="CaloParticles")


# Below this kinetic energy the ratio is dominated by dividing by almost nothing.
REC_RATIO_MIN_EKIN = 0.1      # GeV


def figure_rec_over_ekin(frames, tag):
    """Hit energy (recEnergy) over kinetic energy at the crossing, per SimCluster.
    A shower cannot deposit more than its particle brought across: the legacy
    table peaks just below 1 (calibrated hits), while boundary SimClusters reach
    hundreds -- they hold the hits of folded-in particles but only the root
    track's energy. Each curve is normalised to its own SimClusters, since the
    legacy table has 5-50 times more; kinetic energy, since rest mass is not
    deposited."""
    fig, ax = single_figure()
    bins = np.logspace(-3, 4, 71)
    for key in KEYS:
        d = frames[key][frames[key].ekin > REC_RATIO_MIN_EKIN]
        ratio = (d.rec / d.ekin).to_numpy()
        # SimClusters without hit energy have no place on a log axis: left out of
        # the curve, but counted in its normalisation and in the legend numbers.
        shown = np.clip(ratio[ratio > 0], bins[0], bins[-1] * 0.999)
        counts, edges = np.histogram(shown, bins=bins)
        ax.stairs(counts / max(len(ratio), 1), edges, color=COLORS[key], linestyle=LINESTYLE,
                  label=f"{SERIES[key]}  (median {np.median(ratio):.2f}, "
                        f"{100 * (ratio > 1.5).mean():.0f}% above 1.5)")
    ax.axvline(1.0, ymax=0.62, color=INK_2, linewidth=1.3, linestyle=":")   # stops below the legend
    ax.annotate("hit energy =\ncrossing energy", xy=(1.0, 0.04), xycoords=("data", "axes fraction"),
                xytext=(5, 0), textcoords="offset points", ha="left", va="bottom",
                fontsize=FONT - 4, color=INK_2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    return finish(fig, ax, f"rec_over_ekin_{tag}.png", loc="upper right",
                  title="Hit energy against the energy brought across, per SimCluster",
                  xlabel="Reconstructed energy / $E_\\mathrm{kin}$ at the crossing",
                  ylabel="Fraction of SimClusters")


def figure_rec_energy(rec, e_gun, tag):
    """Summed recEnergy of each CaloParticle's SimClusters. recEnergy is the
    hit energy assigned exclusively to one SimCluster (docs/data_handling.md), so the
    sum over a CaloParticle counts every hit once: the same in both tables."""
    fig, ax = single_figure()
    # The sample's own range, every CaloParticle inside it (0.45-1.74 E_gun on
    # the local samples): no clipping, no empty decade on the left.
    lo = min(rec[k].min() for k in KEYS)
    lo = 0.0 if lo <= 0 else 0.95 * lo
    bins = np.linspace(lo, 1.02 * max(rec[k].max() for k in KEYS), 71)
    for key in KEYS:
        r = rec[key]
        steps(ax, r, bins, key,
              f"{SERIES[key]}  (mean {r.mean():.3g} GeV, median {np.median(r):.3g} GeV)")
    ax.set_yscale("log")
    return finish(fig, ax, f"rec_energy_{tag}.png", loc="upper left",
                  title="Reconstructed energy per CaloParticle",
                  xlabel="$\\Sigma$ reconstructed energy per CaloParticle  [GeV]",
                  ylabel="CaloParticles")


def figure_n_per_caloparticle(cp_counts, tag):
    """One bin per count up to 10, where the boundary table lives (mostly 1-2),
    wider bins beyond for the legacy table's tens to hundreds; the axis is
    linear up to 10 and logarithmic after, and heights are per unit of count so
    the wider bins add no steps."""
    fig, ax = single_figure()
    top = max(int(cp_counts[k].max()) for k in KEYS)
    bins = np.arange(-0.5, 10.5)
    if top > 10:
        bins = np.unique(np.r_[bins, np.round(np.geomspace(10.5, top + 1.5, 25)) + 0.5])
    for key in KEYS:
        n = cp_counts[key]
        counts, edges = np.histogram(n, bins=bins)
        ax.stairs(counts / np.diff(edges), edges, color=COLORS[key], linestyle=LINESTYLE,
                  label=f"{SERIES[key]}  (mean {n.mean():.2f}, median {np.median(n):.0f})")
    if top > 10:
        ax.set_xscale("symlog", linthresh=10, linscale=1.2)
    ticks = [t for t in (0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500) if t <= bins[-1]]
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.set_xticks([], minor=True)
    ax.set_xlim(-0.5, bins[-1])
    ax.set_yscale("log")
    return finish(fig, ax, f"n_per_caloparticle_{tag}.png", title="SimClusters per CaloParticle",
                  xlabel="Number of SimClusters per CaloParticle",
                  ylabel="CaloParticles per unit multiplicity")


def figure_n_per_event(counts, tag):
    fig, ax = single_figure()
    bins = count_bins(max(counts[k].max() for k in KEYS))
    for key in KEYS:
        steps(ax, counts[key], bins, key, f"{SERIES[key]}  (mean {counts[key].mean():.1f})")
    ax.set_yscale("log")
    return finish(fig, ax, f"n_per_event_{tag}.png", title="SimClusters per event",
                  xlabel="Number of SimClusters per event", ylabel="Events")


# --- energy-weighted spectrum -------------------------------------------------------
# Each SimCluster enters with weight E_kin, so a bin's height is the kinetic energy
# carried across by objects of that energy, summed over the sample. `fraction` is
# E_kin / E_gun, the axis the FineCalo excess lives on; `gev` is
# the same bins times E_gun.

# Per table, in units of E_gun.
#   legacy:   the top edge sits past 1 so the last bin is not cut at the
#             primary's energy; the bottom is below the softest crossing either
#             config records at 200 GeV.
#   boundary: cut to where the weight is. Below 3e-3 lies at most 0.33% of the
#             summed E_kin, and nothing in the boundary table exceeds E_gun.
#             What falls outside is printed.
WEIGHTED_BINS = {"legacy": np.logspace(-6, np.log10(1.1), 61),
                 "boundary": np.logspace(np.log10(3e-3), np.log10(1.01), 41)}
# Boundary table: photons put 48-68% of their weight in the bin at E_gun, which on
# a linear axis flattens a tail reaching 0.005 E_gun, so their y axis is log.
LOG_Y = {"photon": True}


def gun_energy(tree):
    # CaloPart_energy is exactly the gun energy in every local file (output_paths).
    return float(np.median(np.concatenate(tree["CaloPart_energy"].array(library="np"))))


def figure_weighted(frames, e_gun, axis, title, tag, bins, logy):
    fig, ax = single_figure()
    scale = 1.0 if axis == "fraction" else e_gun
    bins = bins * scale
    for key in KEYS:
        ekin = frames[key].ekin.to_numpy()
        x = ekin / e_gun * scale
        sums, _ = np.histogram(x, bins=bins, weights=ekin)
        ax.stairs(sums, bins, color=COLORS[key], linestyle=LINESTYLE, label=SERIES[key])
        outside = ekin[(x < bins[0]) | (x >= bins[-1])]
        print(f"  {SERIES[key]}: {len(ekin)} SimClusters, sum E_kin {ekin.sum():.1f} GeV; "
              f"outside the axis: {len(outside)} objects, {outside.sum():.3g} GeV")
    ax.set_xscale("log")
    ax.set_xlim(bins[0], bins[-1])
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(bottom=1.0)
    else:
        ax.set_ylim(bottom=0)
    xlabel = ("$E_\\mathrm{kin}$ at the crossing / $E_\\mathrm{gun}$"
              if axis == "fraction" else "$E_\\mathrm{kin}$ at the crossing  [GeV]")
    return finish(fig, ax, f"energy_weighted_{axis}_{tag}.png", title=title, xlabel=xlabel,
                  ylabel="$\\Sigma$ $E_\\mathrm{kin}$ per bin  [GeV]")


# --- printed summaries ----------------------------------------------------------------


def summarise(frames, counts):
    for key in KEYS:
        d = frames[key]
        plane, mask = on_face(d)
        off = d.abs_z[~mask]
        extent = f"{off.min():.1f}-{off.max():.1f} cm" if len(off) else "none"
        print(f"{LABELS[key]}: {len(d):,} SimClusters over {len(counts[key]):,} events, "
              f"{counts[key].mean():.1f} per event\n"
              f"   front face |z| : {plane:.3f} cm, {100 * mask.mean():.3f}% on it, "
              f"{len(off)} off it ({extent})\n"
              f"   r range        : {d.r.min():.1f} – {d.r.max():.1f} cm\n"
              f"   max |eta|      : {d.eta.abs().max():.2f}\n"
              f"   median E       : {d.energy.median():.4f} GeV\n"
              f"   median pt      : {d.pt.median():.4f} GeV\n"
              f"   photon share   : {100 * (d.pdg == 22).mean():.1f}%")


# --- driver -----------------------------------------------------------------------------

FIGURES = ("position_z", "position_z_off", "eta", "energy", "pt",
           "type_count", "type_momentum", "type_ratio", "n_per_event", "n_per_caloparticle",
           "ekin_per_caloparticle", "rec_energy", "rec_over_ekin", "energy_weighted")


def draw_figures(frames, counts, cp_counts, rec, ekin, only, tag, e_gun, weighted_bins, logy,
                 weighted_title):
    """Every requested figure from `frames` / `counts`, keyed like KEYS; the paths."""
    written = []

    def wrote(path, indent=""):
        print(f"{indent}wrote", path)
        written.append(path)

    simple = {"position_z": figure_position_z, "position_z_off": figure_position_z_off,
              "eta": figure_eta, "energy": figure_energy, "pt": figure_pt,
              "type_count": figure_type_count, "type_momentum": figure_type_momentum,
              "type_ratio": figure_type_ratio, "rec_over_ekin": figure_rec_over_ekin}
    for name in only:
        if name == "n_per_event":
            wrote(figure_n_per_event(counts, tag))
        elif name == "n_per_caloparticle":
            wrote(figure_n_per_caloparticle(cp_counts, tag))
        elif name == "rec_energy":
            wrote(figure_rec_energy(rec, e_gun, tag))
        elif name == "ekin_per_caloparticle":
            wrote(figure_ekin_per_caloparticle(ekin, e_gun, tag))
        elif name == "energy_weighted":
            for key in KEYS:
                n = counts[key]
                print(f"  {SERIES[key]}: {n.mean():.2f} per event (max {n.max()}), "
                      f"sum E_kin {frames[key].ekin.sum() / len(n):.2f} GeV per event "
                      f"(2 x E_gun = {2 * e_gun:.0f})")
            for axis in ("fraction", "gev"):
                print(f"energy_weighted {axis}:")
                wrote(figure_weighted(frames, e_gun, axis, weighted_title, tag, weighted_bins, logy),
                      indent="  ")
        else:
            wrote(simple[name](frames, tag))
    return written


def weighted_title(sample):
    return (f"Energy-weighted boundary spectrum, {PARTICLE_TEX.get(sample.particle, sample.particle)} "
            f"{sample.energy.replace('GeV', ' GeV')}")


def run(sample, table=None, only=FIGURES, tag=None):
    """Every requested figure for one sample; prints the numbers it drew and
    returns the paths of the figures it wrote."""
    global FIGDIR, P, COLORS, SERIES, LINESTYLE, TABLE_NOW, ENERGY_NOW
    table = table or sample.table
    if table not in sample.tables:
        raise SystemExit(f"{table} table not in both files ({sample.tables})")
    COLORS = config_colours(table)
    SERIES = {k: series_label(k, table, sample.energy) for k in KEYS}
    LINESTYLE = energy_linestyle(sample.energy)
    TABLE_NOW, ENERGY_NOW = table, sample.energy.replace("GeV", " GeV")
    P = samples.TABLES[table]
    tag = (tag or sample.tag) + ("_bound" if table == "boundary" else "")
    trees = {k: sample.tree(k) for k in KEYS}
    FIGDIR = table_dir(trees["fc0"], table)

    print(f"sample {sample.tag}: {sample.particle} {sample.energy}, {table} table")
    for key in samples.KEYS:
        print(f"  {LABELS[key]:9s} {sample.files[key].name}")
    print()

    use_house_style()
    frames = {k: samples.simcluster_frame(trees[k], table, k) for k in KEYS}
    counts = {k: per_event(frames[k], trees[k].num_entries) for k in KEYS}
    cp_counts = {k: per_caloparticle(frames[k], trees[k]) for k in KEYS}
    rec = {k: per_caloparticle(frames[k], trees[k], "rec") for k in KEYS}
    ekin = {k: per_caloparticle(frames[k], trees[k], "ekin") for k in KEYS}
    summarise(frames, counts)
    print()
    logy = table == "boundary" and LOG_Y.get(sample.particle, False)
    return draw_figures(frames, counts, cp_counts, rec, ekin, only, tag, gun_energy(trees["fc0"]), WEIGHTED_BINS[table],
                        logy, weighted_title(sample))


def compare(sample, only=FIGURES, tag=None):
    """The same figures with both SimCluster collections overlaid -- four series,
    one colour each (`hep_style.CONFIG_COLLECTION`) -- into
    outputs/<particle>/<energy>/comparisons/<name>_<tag>_compare.png. Needs a
    pair carrying both tables; returns the paths."""
    global KEYS, FIGDIR, P, COLORS, SERIES, LINESTYLE, ENERGY_NOW
    if not {"boundary", "legacy"} <= set(sample.tables):
        raise SystemExit(f"{sample.tag} does not carry both SimCluster tables ({sample.tables})")
    saved = KEYS, FIGDIR, P, COLORS, SERIES, LINESTYLE, ENERGY_NOW
    try:
        KEYS = COMPARE_KEYS
        ENERGY_NOW = sample.energy.replace("GeV", " GeV")
        COLORS = {k: CONFIG_COLLECTION[k] for k in KEYS}
        SERIES = {k: series_label(k[0], k[1], sample.energy) for k in KEYS}
        LINESTYLE = energy_linestyle(sample.energy)
        P = "SimCluster"                  # both tables: the branch names differ in prefix only
        trees = {k: sample.tree(k) for k in samples.KEYS}
        FIGDIR = style_dir(trees["fc0"], "comparisons")
        tag = (tag or sample.tag) + "_compare"
        print(f"sample {sample.tag}: {sample.particle} {sample.energy}, both tables -> {FIGDIR}")
        use_house_style()
        frames = {k: samples.simcluster_frame(trees[k[0]], k[1], k[0]) for k in KEYS}
        counts = {k: per_event(frames[k], trees[k[0]].num_entries) for k in KEYS}
        cp_counts = {k: per_caloparticle(frames[k], trees[k[0]]) for k in KEYS}
        rec = {k: per_caloparticle(frames[k], trees[k[0]], "rec") for k in KEYS}
        ekin = {k: per_caloparticle(frames[k], trees[k[0]], "ekin") for k in KEYS}
        # Legacy bins span both tables; a log y axis keeps the boundary table's
        # few, energetic objects and the legacy table's many soft ones in view.
        return draw_figures(frames, counts, cp_counts, rec, ekin, only, tag, gun_energy(trees["fc0"]),
                            WEIGHTED_BINS["legacy"], True, weighted_title(sample))
    finally:
        KEYS, FIGDIR, P, COLORS, SERIES, LINESTYLE, ENERGY_NOW = saved
