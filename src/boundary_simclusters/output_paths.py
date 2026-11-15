"""Where a figure belongs in `outputs/`, derived from the data it came from.

`outputs/` is organised by **particle type, then gun energy, then SimCluster
table**:

    outputs/<particle>/<energy>/<table>/<name>_<tag>.png     # table_dir
    outputs/<particle>/<energy>/<family>/<name>_<tag>.png    # style_dir, no table
    outputs/slides/<table>/<name>.png                        # slide_path, cross-sample

so `outputs/photon/20GeV/boundary/n_per_event_g20_bound.png` next to
`outputs/photon/20GeV/legacy/n_per_event_gamma20.png`. The table level is the
collection the figure was drawn from — `boundary` (`BoundarySimCluster_*`) or
`legacy` (`SimCluster_*`) — so the two collections of one sample never share a
folder. The standard run and `event_scatter.py` write to
the table folder. Figures tied to no table use a family folder instead:
`envelope` (`classify_boundary_surface.py`), `comparisons`, `collections`.

The sample tag stays in the filename even though the path repeats it. It is
redundant on disk and load-bearing off it: figures get dragged into slide decks,
where three files called `energy.png` are indistinguishable.

**Particle and energy are read from the ROOT file, never parsed from its name.**
The filename prefixes in this repo are production labels that do not match the
pdgId — `11_` is the photon gun and `22_` the electron gun (see README). A path
built from the filename would file half the samples under the wrong particle.
"""

from contextlib import contextmanager
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
# Where everything is written. Code reads these as `output_paths.OUT` /
# `output_paths.PROCESSED` when it writes, never copies them at import, so the
# tests can point both at a temporary folder.
OUT = REPO / "outputs"
PROCESSED = REPO / "data" / "processed"

# Only what this project's guns actually fire. An unknown code is labelled
# rather than guessed at, so a new sample lands somewhere obvious and wrong
# instead of somewhere plausible and wrong.
PARTICLE_NAMES = {
    11: "electron", -11: "positron",
    22: "photon",
    13: "muon", -13: "muon",
    211: "pion", -211: "pion",
    2212: "proton", 2112: "neutron",
}


def particle_name(tree):
    """Gun particle, from GenPart_pdgId. Mixed guns get a joined name."""
    pdg = np.concatenate(tree["GenPart_pdgId"].array(library="np"))
    codes = sorted({int(c) for c in pdg})
    names = []
    for code in codes:
        name = PARTICLE_NAMES.get(code, f"pdg{code}")
        if name not in names:
            names.append(name)
    return "_".join(names)


def energy_label(tree):
    """Gun energy as a directory name, e.g. "200GeV".

    `CaloPart_energy` is exactly the gun energy in every local file (200.000 /
    20.000), so its median is the sample's nominal energy. A non-integer energy
    is written with `p` for the point, since `.` in a directory name reads as an
    extension.
    """
    energy = float(np.median(np.concatenate(
        tree["CaloPart_energy"].array(library="np"))))
    if abs(energy - round(energy)) < 0.01 * max(energy, 1.0):
        return f"{round(energy)}GeV"
    return f"{energy:.1f}".replace(".", "p") + "GeV"


def sample_dir(tree):
    """`outputs/<particle>/<energy>` for the sample this tree holds."""
    return OUT / particle_name(tree) / energy_label(tree)


def style_dir(tree, style, table=None):
    """`outputs/<particle>/<energy>/[<table>/]<style>`, created if missing.

    Scripts resolve this once from their fc0 tree and save every figure
    against it, so the routing decision is made in one place per run.
    `table` is "boundary" or "legacy", the SimCluster collection drawn.
    """
    directory = sample_dir(tree) / (table or "") / style
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def table_dir(tree, table):
    """`outputs/<particle>/<energy>/<table>`, created if missing: the one folder
    the standard run writes a sample's figures for one SimCluster table into."""
    directory = sample_dir(tree) / table
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def slide_path(name, table=None):
    """Full path for a cross-sample presentation figure, under
    `outputs/slides/[<table>/]`."""
    directory = OUT / "slides" / (table or "")
    directory.mkdir(parents=True, exist_ok=True)
    return directory / name


@contextmanager
def _writing_to(out, processed):
    global OUT, PROCESSED
    saved = OUT, PROCESSED
    OUT, PROCESSED = Path(out), Path(processed)
    try:
        yield
    finally:
        OUT, PROCESSED = saved


def redirected(workdir):
    """Point `OUT` and `PROCESSED` at `workdir/outputs`, `workdir/processed` for
    the duration, so a run writes nothing into the repo (the tests)."""
    return _writing_to(Path(workdir) / "outputs", Path(workdir) / "processed")


def rooted(root):
    """Point `OUT` and `PROCESSED` at `root/outputs`, `root/data/processed` for
    the duration: an analysis tree laid out like the repo's own, for samples kept
    apart from the main ones (a second production, say `other/`). Keep such files
    in their own folder: partner lookup goes by GEN overlap, so two productions of
    the same gun in one folder can pair across productions."""
    return _writing_to(Path(root) / "outputs", Path(root) / "data" / "processed")
