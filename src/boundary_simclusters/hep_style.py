"""House figure style: light typography in a dark frame, one panel per figure.

The lighter look of the slide figures -- small type, muted ink, wide panels --
carrying the four fixes worked out against the ROOT-style draft:

- **Legends never overlap the data.** `overlay()` measures the rendered legend
  and clears exactly that much axis above the tallest bin, so no series can run
  underneath it, whatever the labels say.
- **A dark, closed frame**, rather than the slide style's pale open one.
- **Ticks, not lines.** Tick marks point inward on all four sides and there are
  no grid lines: numbers are carried by the frame, never by lines drawn across
  the data.
- **Axis labels centred**, not ROOT's right/top-anchored convention.

Colours match `pipeline.make_figures`, so one colour means one
config everywhere in the repo. Only the font (Tex Gyre Heros) is still inherited from mplhep.

Figures name their series `FineCalo` and `Standard` -- never `fc0`/`fc1`, which
mean nothing outside this repo.

`save()` writes every figure as a PNG, as drawn: titles, no experiment label.
A figure from these samples is not CMS-approved. If you add a label, it says
"Private work (CMS simulation)", never "CMS Preliminary", nor "CMS" or
"CMS Simulation" alone. The PDF version is switched off in `save()`.

Usage:
    from boundary_simclusters.hep_style import use_house_style, single_figure, overlay, save
    use_house_style()
    fig, ax = single_figure()
    overlay(ax, {"Standard": a, "FineCalo": b}, bins, xlabel="energy [GeV]",
            colors=config_colours("boundary"))
    save(fig, path)

The ratio panel (`ratio_figure`/`compare`) was removed 2026-10-11, unused: on a
log axis the gap between two curves already is their ratio. It is in the git
history if a linearly plotted variable ever needs it.
"""

# from pathlib import Path    # for the PDF version, switched off in save()

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

# Only the generic PALETTE below uses this pair. Figures take their config
# colours from CONFIG_COLLECTION.
STANDARD, FINECALO = "#2a78d6", "#eb6834"

# Colour identifies config *and* SimCluster collection (decided 2026-10-11):
# purple / gold for Standard / FineCalo on the boundary collection, blue / red
# on the legacy one -- in figures of one collection (`config_colours(table)`)
# and in figures crossing both (four series). Validated
# with the data-viz palette checks, every pair and as adjacent lines: worst
# normal-vision Delta E 16.3 (floor 15), worst colour-blind 13.0 (target 8).
# Gold is light on white (2.1:1), so these figures keep the legend and the line
# style (solid boundary, dashed legacy) as a second cue.
CONFIG_COLLECTION = {("fc0", "boundary"): "#4a3aa7", ("fc0", "legacy"): "#2a78d6",
                     ("fc1", "boundary"): "#eda100", ("fc1", "legacy"): "#e34948"}


def config_colours(table):
    """{"fc0": colour, "fc1": colour} for a figure of one SimCluster collection."""
    return {key: CONFIG_COLLECTION[(key, table)] for key in ("fc0", "fc1")}


# Figure rules (decided 2026-10-11): every legend entry names the SimCluster
# collection and the gun energy, even in a figure of one collection or one
# energy; axis labels say
# "SimClusters", never "boundary SimClusters"; the line style is the gun
# energy -- solid 200 GeV, dashed 20 GeV -- in every figure with lines.
CONFIG_NAMES = {"fc0": "Standard", "fc1": "FineCalo"}
ENERGY_LINESTYLE = {"200GeV": "-", "20GeV": "--"}


def energy_text(energy):
    """ "200GeV" -> "200 GeV"; a list of labels -> "20 and 200 GeV" (ascending)."""
    if isinstance(energy, str):
        return energy.replace("GeV", " GeV")
    values = sorted({e.replace("GeV", "") for e in energy}, key=lambda v: float(v.replace("p", ".")))
    return (", ".join(values[:-1]) + " and " + values[-1] if len(values) > 1 else values[0]) + " GeV"


def series_label(key, table, energy):
    """Legend name of one config's series: "Standard, boundary, 200 GeV". `energy`
    is a label ("200GeV") or, for a series spanning samples, a list of them."""
    return f"{CONFIG_NAMES[key]}, {table}, {energy_text(energy)}"


def energy_linestyle(energy):
    """Line style of a sample by its energy label ("200GeV", "20GeV"); solid otherwise."""
    return ENERGY_LINESTYLE.get(energy, "-")

# Species, for figures whose series are particle types rather than configs
# (the topology slides): green electrons, magenta photons -- hues the config
# schemes above do not use.
ELECTRON_COLOUR, PHOTON_COLOUR = "#008300", "#e87ba4"
PALETTE = (FINECALO, STANDARD, "#e42536", "#964a8b", "#7a21dd", "#9c9ca1")

# Ink scale from the slide figures. All text is set in INK: the grey steps
# below it are kept for callers that want a recessive annotation, but grey type
# on a small font was hard to read on a projector and in print.
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#898781"
FRAME = INK
GUIDE = "#b4b2ab"

# Colour per particle species, for the figures that encode species rather than
# configuration. Deliberately disjoint from the config colours, which must
# not also come to mean a particle. Labels are mathtext because the house font has no U+207A/U+207B.
SPECIES = {
    22: (r"$\gamma$", "#7a21dd"),
    11: (r"$e^-$", "#e42536"),
    -11: (r"$e^+$", "#2e9e5b"),
    13: (r"$\mu^-$", "#964a8b"),
    -13: (r"$\mu^+$", "#964a8b"),
}
SPECIES_OTHER = ("other", "#0f8b8d")


def species(pdg):
    """(mathtext label, colour) for a pdgId. Anything else is 'other'.

    Hadrons collapse into one entry on purpose: in every boundary population
    measured here they are together well under 0.2%, so
    naming each would spend a legend row on a species with no markers in it.
    """
    return SPECIES.get(int(pdg), SPECIES_OTHER)

# Large, so figures stay readable when printed at half the text width.
BASE_FONT = 24

# Fraction of the axes height left free above the tallest bin, so the legend
# has somewhere to sit that no series can reach.
LEGEND_HEADROOM = 0.16
# Smallest legend type before an entry is split over two lines.
LEGEND_MIN_FONT = 14


def use_house_style():
    """Install the house look globally. Call once, before creating any figure."""
    # mplhep is kept only for its font; everything visual below is overridden.
    plt.style.use(hep.style.ROOT)
    plt.rcParams.update({
        "figure.dpi": 200,
        "savefig.dpi": 200,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.prop_cycle": plt.cycler(color=PALETTE),
        "axes.grid": False,
        "axes.linewidth": 1.4,
        "axes.edgecolor": FRAME,
        "axes.titlesize": BASE_FONT + 2,
        "axes.titlecolor": INK,
        "axes.labelsize": BASE_FONT,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.labelsize": BASE_FONT - 2,
        "ytick.labelsize": BASE_FONT - 2,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "legend.fontsize": BASE_FONT - 4,
        "legend.frameon": False,
        "legend.labelcolor": INK,
        "lines.linewidth": 1.6,
        # Ticks inward on all four sides; no grid anywhere.
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 5, "ytick.major.size": 5,
        "xtick.minor.size": 3, "ytick.minor.size": 3,
        "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "xtick.minor.width": 0.8, "ytick.minor.width": 0.8,
        # mplhep's ROOT style right-aligns x labels and top-anchors y labels.
        "xaxis.labellocation": "center", "yaxis.labellocation": "center",
    })


def _frame(ax):
    for spine in ax.spines.values():
        spine.set_color(FRAME)
    return ax


def single_figure(figsize=(9.0, 6.4)):
    """One framed panel, for distributions that have no meaningful ratio."""
    fig, ax = plt.subplots(figsize=figsize)
    return fig, _frame(ax)


def panel_figure(ncols=2, figsize=(13.0, 6.4)):
    """A row of independently framed panels sharing no axes."""
    fig, axes = plt.subplots(1, ncols, figsize=figsize)
    for ax in np.atleast_1d(axes):
        _frame(ax)
    return fig, axes


def _reserve_headroom(ax, fraction=LEGEND_HEADROOM):
    """Grow the y-axis so the top `fraction` of the panel holds no data."""
    fraction = min(max(fraction, 0.0), 0.38)
    lo, hi = ax.get_ylim()
    if ax.get_yscale() == "log":
        ax.set_ylim(lo, lo * 10 ** (np.log10(hi / lo) / (1 - fraction)))
    else:
        ax.set_ylim(lo, lo + (hi - lo) / (1 - fraction))


def legend_with_headroom(ax, loc="upper right", pad=0.02, **kw):
    """Place the legend, then clear exactly as much axis as it actually needs.

    Measured rather than assumed: a multi-line legend is taller than any fixed
    reserve, and guessing produced overlaps on precisely the plots whose labels
    carry the most information. The legend is anchored in axes coordinates, so
    growing the y-axis afterwards moves the data down and leaves the legend
    where it is -- one pass is enough.
    """
    legend = ax.legend(loc=loc, **kw)
    try:
        renderer = ax.figure.canvas.get_renderer()
    except AttributeError:                      # backend without a live renderer
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
    box = legend.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    # Wider than the panel: shrink the
    # legend type a point at a time so every entry stays on one line; only
    # below LEGEND_MIN_FONT move the "(mean ...)" part onto a second line.
    size = float(kw.pop("fontsize", plt.rcParams["legend.fontsize"]))
    while (box.x0 < 0 or box.x1 > 1) and size > LEGEND_MIN_FONT:
        size -= 1
        legend.remove()
        legend = ax.legend(loc=loc, fontsize=size, **kw)
        box = legend.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    if box.x0 < 0 or box.x1 > 1:
        handles, labels = ax.get_legend_handles_labels()
        legend.remove()
        legend = ax.legend(handles, [t.replace("  (", "\n(", 1) for t in labels], loc=loc,
                           fontsize=size, **kw)
        box = legend.get_window_extent(renderer).transformed(ax.transAxes.inverted())
    _reserve_headroom(ax, max(LEGEND_HEADROOM, 1.0 - box.y0 + pad))
    return legend


def overlay(ax, series, bins, *, title=None, xlabel=None, ylabel="Counts",
            logy=True, colors=None, linestyles=None, legend_loc="upper right",
            legend_kw=None):
    """Step histograms of every entry in `series` (label -> values) on one axes.

    Returns label -> bin counts, so callers can quote numbers without
    re-histogramming.

    `linestyles` is for figures carrying a second dimension beyond the config --
    an energy scan overlays four series, and colour is already spoken for
    (config and collection, `CONFIG_COLLECTION`). Passing a dashed style for the
    second sample keeps colour meaning what it means everywhere else.
    """
    colors = colors or PALETTE
    linestyles = linestyles or ["-"] * len(series)
    counts = {}
    for (label, values), colour, dash in zip(series.items(), colors, linestyles):
        counts[label], _ = np.histogram(np.asarray(values), bins=bins)
        hep.histplot(counts[label], bins, ax=ax, histtype="step",
                     color=colour, linestyle=dash, label=label)

    if logy:
        ax.set_yscale("log")
        # Empty bins would otherwise pull the axis down to ~1e-1 and waste a
        # decade below the lowest real count.
        ax.set_ylim(bottom=0.5)
    ax.set_xlim(float(np.asarray(bins)[0]), float(np.asarray(bins)[-1]))
    ax.set_ylabel(ylabel, linespacing=1.5)
    if xlabel is not None:
        ax.set_xlabel(xlabel, linespacing=1.5)
    if len(series) > 1:
        legend_with_headroom(ax, loc=legend_loc, **(legend_kw or {}))
    if title:
        ax.set_title(title, pad=12)
    return counts


# Only the PDF version in save(), switched off, calls this and the label code below.
def _pdf_version(fig):
    """Turn a drawn figure into its PDF form: no titles, the private-work label.

    The figure carries no title -- a caption does that job -- and says
    "Private work (CMS simulation)", never
    "CMS", "CMS Preliminary" or "CMS Simulation" alone. The label sits inside
    the top-left corner of the first panel's frame; `_private_work_label` keeps
    data and legend clear of it. Returns the titles removed, joined by " | ",
    for the PDF's metadata.
    """
    titles = []
    for ax in fig.axes:
        for loc in ("center", "left", "right"):
            if ax.get_title(loc=loc):
                titles.append(" ".join(ax.get_title(loc=loc).split()))
                ax.set_title("", loc=loc)
    if fig.axes:
        _private_work_label(fig.axes[0])
    return " | ".join(titles)


PRIVATE_WORK = "Private work (CMS simulation)"


def _private_work_label(ax, pad=0.02):
    """ "Private work (CMS simulation)", all italic, inside the top-left corner of
    `ax`'s frame.

    The figure is already drawn, so the label has to make room for itself: a
    legend it would touch moves down below it, and the y-axis grows (as in
    `legend_with_headroom`) until no data reaches the label or the legend.
    """
    label = ax.text(0.025, 0.975, PRIVATE_WORK, transform=ax.transAxes, ha="left", va="top",
                    fontsize=BASE_FONT - 2, fontstyle="italic", color=INK)
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_axes = ax.transAxes.inverted()
    box = label.get_window_extent(renderer).transformed(to_axes)
    top = box.y0 - pad                               # highest point data may reach
    legend = ax.get_legend()
    if legend is not None:
        leg = legend.get_window_extent(renderer).transformed(to_axes)
        if leg.x0 < box.x1 and leg.y1 > box.y0:      # the two would overlap
            legend.set_bbox_to_anchor((0.0, 0.0, 1.0, top), transform=ax.transAxes)
            fig.canvas.draw()
            leg = legend.get_window_extent(renderer).transformed(to_axes)
        if leg.y1 > 0.75:                            # a legend at the top caps the data
            top = min(top, leg.y0 - pad)
    if top > 0.2:
        _clear_below(ax, top)
        _clear_texts_below(ax, top, skip=label)


def _clear_texts_below(ax, top, skip=None, rounds=6):
    """Like `_clear_below`, for text placed in data coordinates (value labels
    above bars, say), which the data limits do not include: grow the y-axis
    until the highest such text ends below `top` (axes fraction)."""
    to_axes = ax.transAxes.inverted()
    for _ in range(rounds):
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
        tops = [t.get_window_extent(renderer).transformed(to_axes).y1 for t in ax.texts
                if t is not skip and t.get_visible() and t.get_text()
                and getattr(t, "xycoords", None) in ("data", None)
                and (getattr(t, "xycoords", None) == "data" or t.get_transform() == ax.transData)]
        highest = max(tops, default=0.0)
        if highest <= top:
            return
        lo, hi = ax.get_ylim()
        grow = highest / top * 1.02
        if ax.get_yscale() == "log" and lo > 0:
            ax.set_ylim(lo, lo * (hi / lo) ** grow)
        else:
            ax.set_ylim(lo, lo + (hi - lo) * grow)


def _clear_below(ax, top):
    """Grow the y-axis just enough that the highest data point sits at or below
    `top` (axes fraction). Measured on the data limits, so headroom a legend has
    already reserved is counted rather than added to."""
    lo, hi = ax.get_ylim()
    y1 = ax.dataLim.y1
    if not np.isfinite(y1) or y1 <= lo:
        return
    if ax.get_yscale() == "log":
        if lo <= 0:
            return
        need = lo * 10 ** (np.log10(y1 / lo) / top)
    else:
        need = lo + (y1 - lo) / top
    if need > hi:
        ax.set_ylim(lo, need)


def save(fig, path):
    """Write `path` (the PNG, the figure as drawn: titles, no experiment label); return `path`."""
    fig.savefig(path, bbox_inches="tight")
    # The PDF version, with no title and the private-work label, is switched off.
    # To write it, uncomment these lines and the pathlib import at the top.
    # title = _pdf_version(fig)
    # fig.savefig(Path(path).with_suffix(".pdf"), bbox_inches="tight",
    #             metadata={"Title": title} if title else None)
    plt.close(fig)
    return path

