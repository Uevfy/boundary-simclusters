"""How far each config's boundary sits from the calorimeter, and
what Standard's interior SimClusters are. A standalone study, not part of
the standard run.

The face classification this script used to do -- against surface radii
calibrated into `data/processed/boundary_surfaces.json` -- now lives in
the surfaces step of `analyze.py`, against the CMSSW source volumes; the two agree on every
SimCluster but 2-3 per FineCalo sample, at the corner where the front face meets
the inner bore. What stays here is what a census cannot show:

1. **The envelope, measured.** The front face and every cylindrical surface,
   calibrated from this pair's SimTrack crossings (`calibrate_surfaces`) and
   printed beside the source radii, with the containment and staircase checks
   that show the faces form one closed shape.
2. **The standoff** -- how far each face sits from the nearest sensitive
   material, using RecHits as the map of where the calorimeter actually is.
   This is what shows Standard recording 28.65 cm upstream of anything
   calorimetric while FineCalo records flush against it.
3. **Interior SimClusters**: per CaloParticle whose every SimCluster is interior,
   did the primary cross the front face, and did that crossing get a SimCluster
   of its own? This is Standard recording tracks born inside CALO at their
   vertex.

Reads the legacy table (`SimCluster_*`);
no event pairing is needed.

Usage:
    python scripts/classify_boundary_surface.py testNanoML.root FinetestNanoML.root
    python scripts/classify_boundary_surface.py 11_fc1_s1_testNanoML.root --tag gamma_s1
"""

import argparse
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from boundary_simclusters import samples  # noqa: E402
from boundary_simclusters.boundary_geometry import (ENVELOPE, SURFACE_TOL, calibrate_surfaces,  # noqa: E402
                               classify_faces, front_face, front_face_mask,
                               radius, surface_name)
from boundary_simclusters.hep_style import (INK, INK_2, config_colours, energy_linestyle,  # noqa: E402
                                            series_label,
                       panel_figure, save, use_house_style)
from boundary_simclusters.output_paths import style_dir  # noqa: E402


KEYS = ("fc0", "fc1")
LABELS = {"fc0": "Standard", "fc1": "FineCalo"}
COLORS = config_colours("legacy")         # this study reads the legacy table

# RecHits are identical in both configs (docs/data_handling.md), so the material map is
# read from fc0 only, and from a slice of it: 300 events is already 13.2M hits,
# and the map is a geometric envelope that does not sharpen with more.
MATERIAL_EVENTS = 300

def load(tree):
    """(|z|, r) of every SimCluster crossing of one config."""
    def flat(branch):
        return np.concatenate(tree[branch].array(library="np"))

    x, y, z = (flat(f"SimCluster_impactPoint_{a}") for a in "xyz")
    return {"abs_z": np.abs(z.astype(np.float64)), "r": radius(x, y)}


def crossing_tracks(tree):
    """Boundary crossings at SimTrack level: (|z|, r) of every track that crossed.

    SimTracks are what the surfaces are *calibrated* from, even though
    SimClusters are what gets classified. A SimCluster additionally requires an
    HGCAL hit, so the SimTrack collection holds 361,550 crossings
    against 160,940 — and 34x as many off the front face, where the surfaces
    live. Both see identical radii wherever both see a surface at all; SimTrack
    simply reaches the faint ones. The r = 38.36 cm bore step is invisible in the
    SimCluster sample of some pairs and unmissable here.
    """
    def flat(branch):
        return np.concatenate(tree[branch].array(library="np"))

    crossed = flat("SimTrack_crossedBoundary").astype(bool)
    x, y, z = (flat(f"SimTrack_boundaryPos_{a}")[crossed] for a in "xyz")
    return np.abs(z.astype(np.float64)), radius(x, y)


def material_map(tree, n_events):
    """Where the sensitive calorimeter material is, from RecHit positions.

    RecHits mark *sensitive* layers only — HGCAL's absorber, support and cooling
    leave none — so this is the front of the instrumented volume, not the front
    of the physical detector. That distinction is why the FineCalo standoff below
    is reported as "flush", not as a measured air gap.
    """
    def flat(branch):
        return np.concatenate(
            tree[branch].array(library="np", entry_stop=n_events))

    x, y, z = (flat(f"RecHitHGC_{a}") for a in "xyz")
    return {"abs_z": np.abs(z.astype(np.float64)), "r": radius(x, y),
            "n_events": n_events}


def section(title):
    print(f"\n{title}\n{'-' * len(title)}")


def source_radii(config):
    """Radii of the source volume's cylindrical faces (constant r along z)."""
    return sorted({round(r0, 2) for _name, z0, r0, z1, r1, *_ in ENVELOPE[config]
                   if z0 != z1 and r0 == r1})


def calibrate(data, tracks):
    """Surface radii per config, measured from this pair's SimTrack crossings.

    SimTracks rather than SimClusters: a SimCluster also needs an HGCAL hit,
    which removes most off-face crossings -- the r = 38.36 cm bore step
    is invisible in the SimCluster sample of some pairs and unmissable here.
    The source radii are printed beside them, so a face this pair is too small
    to reach shows up as missing rather than silently.
    """
    section("1. Boundary surfaces")
    surfaces = {}
    for key in KEYS:
        track_z, track_r = tracks[key]
        off = ~front_face_mask(track_z, front_face(track_z))
        surfaces[key] = calibrate_surfaces(track_r[off])
        d = data[key]
        plane = front_face(d["abs_z"])
        print(f"{LABELS[key]:9s}: front face |z| = {plane:.3f} cm")
        print(f"           surveyed from {int(off.sum()):,} off-face SimTrack crossings "
              f"(SimCluster would give {int((~front_face_mask(d['abs_z'], plane)).sum()):,})")
        print(f"           measured surfaces: {[round(v, 3) for v in surfaces[key]] or 'none'}")
        print(f"           source cylinders:  {source_radii(key)}")
    return surfaces


def report_envelope(tracks, surfaces):
    """The shape: front-face extent, each surface, and the corners that join them.

    Surveyed from SimTracks, not SimClusters. This section describes *geometry*,
    and the SimCluster sample is too thin to measure it: its front-face outer
    edge lands 1.2 cm short of the true one on the 20 GeV pair and 10.2 cm short
    on the photon pair, purely for want of a crossing out there. The SimCluster
    counts per surface are section 4's job.
    """
    section("2. The envelope")
    for key in KEYS:
        track_z, track_r = tracks[key]
        plane = front_face(track_z)
        on_face = front_face_mask(track_z, plane)
        face_r = track_r[on_face]
        d = {"abs_z": track_z, "r": track_r}
        print(f"\n{LABELS[key]}  ({len(track_z):,} SimTrack crossings)")
        print(f"  front face : |z| = {plane:.3f} cm, an annulus "
              f"r = {face_r.min():.3f} to {face_r.max():.2f} cm")
        if not surfaces[key]:
            print("  no other surface reached by this sample")
            continue
        rows = []
        for surface_r in surfaces[key]:
            on = ~on_face & (np.abs(d["r"] - surface_r) < SURFACE_TOL)
            if not on.any():
                rows.append({"surface": surface_name(surface_r, np.median(face_r)),
                             "crossings": 0, "r mean": round(surface_r, 3),
                             "r std": np.nan, "|z| from": np.nan, "|z| to": np.nan})
                continue
            rows.append({
                "surface": surface_name(surface_r, float(np.median(face_r))),
                "crossings": int(on.sum()),
                "r mean": round(float(d["r"][on].mean()), 3),
                "r std": float(f"{d['r'][on].std():.1e}"),
                "|z| from": round(float(d["abs_z"][on].min()), 1),
                "|z| to": round(float(d["abs_z"][on].max()), 1),
            })
        print(pd.DataFrame(rows).set_index("surface").to_string())

        # Containment, not edge-matching: does this surface bound the front face?
        # `face_r.max()` is set by a single stray crossing -- 209 cm on the photon
        # pair against a real edge at 123.3 -- so the test is what share of the
        # face lies *past* the surface, which a handful of outliers cannot move.
        # This is the evidence the faces form one closed shape.
        median_r = float(np.median(face_r))
        inner = [s for s in surfaces[key] if s < median_r]
        outer = [s for s in surfaces[key] if s >= median_r]
        for adjacent, which in ((min(inner, default=None), "inner"),
                                (max(outer, default=None), "outer")):
            if adjacent is None:
                continue
            beyond = (face_r < adjacent) if which == "inner" else (face_r > adjacent)
            share = 100.0 * float(np.mean(beyond))
            verdict = "bounds the front face" if share < 0.1 else "does NOT bound it"
            print(f"  containment: {which} surface r = {adjacent:8.3f} — "
                  f"{share:.4f}% of front-face crossings lie beyond it, {verdict}")
        for below, above in zip(inner, inner[1:]):
            step_z = d["abs_z"][~on_face & (np.abs(d["r"] - above) < SURFACE_TOL)]
            prev_z = d["abs_z"][~on_face & (np.abs(d["r"] - below) < SURFACE_TOL)]
            if not step_z.size or not prev_z.size:
                continue
            print(f"  staircase: r {below:.3f} -> {above:.3f} cm "
                  f"(+{above - below:.3f}) between |z| {prev_z.max():.1f} "
                  f"and {step_z.min():.1f} cm")


def report_standoff(data, tracks, surfaces, material):
    """How far each face sits from the nearest sensitive material."""
    section("3. Standoff from the calorimeter")
    mz, mr = material["abs_z"], material["r"]
    print(f"sensitive material ({material['n_events']} events, {len(mz):,} RecHits): "
          f"|z| = {mz.min():.2f} to {mz.max():.2f} cm, "
          f"r = {mr.min():.2f} to {mr.max():.2f} cm")

    rows = []
    for key in KEYS:
        track_z, track_r = tracks[key]
        d = {"abs_z": track_z, "r": track_r}
        plane = front_face(track_z)
        on_face = front_face_mask(track_z, plane)
        rows.append({
            "config": LABELS[key], "surface": "front face",
            "at": f"|z| = {plane:.2f}",
            "nearest material": f"|z| = {mz.min():.2f}",
            "standoff (cm)": round(float(mz.min() - plane), 2),
            "sense": "upstream of the material" if plane < mz.min() else "inside",
        })
        for surface_r in surfaces[key]:
            on = ~on_face & (np.abs(d["r"] - surface_r) < SURFACE_TOL)
            if not on.any():
                continue
            # Only material at the depths this surface actually spans: the bore
            # steps outward with z, so the global r-min would understate every
            # step but the first.
            span = (mz >= d["abs_z"][on].min()) & (mz <= d["abs_z"][on].max())
            if not span.any():
                continue
            inner = surface_r < np.median(d["r"][on_face])
            edge = mr[span].min() if inner else mr[span].max()
            rows.append({
                "config": LABELS[key],
                "surface": surface_name(surface_r, float(np.median(d["r"][on_face]))),
                "at": f"r = {surface_r:.2f}",
                "nearest material": f"r = {edge:.2f}",
                "standoff (cm)": round(float(edge - surface_r if inner
                                             else surface_r - edge), 2),
                "sense": "inward of the material" if inner else "outward of it",
            })
    print(pd.DataFrame(rows).set_index(["config", "surface"]).to_string())
    print("\nRecHits mark sensitive layers only — absorber and support leave none —")
    print("so a standoff of ~1 cm means flush with the entrance, not a measured gap.")


def report_interior(trees):
    """Crossings strictly inside the envelope, and what their ancestor did.

    An `interior` crossing lies on no face, so no surface was crossed there. The
    Standard mechanism flags a track whose *first step* is already inside CALO,
    with its production vertex as the boundary position. This
    section asks whether that is what the interior population is.

    The schema carries no step number or production vertex, so the test is
    indirect, per CaloParticle: if *every* SimCluster of a CaloParticle is
    interior, did a track carrying most of the primary's energy cross the front
    face on that side, and did that crossing get a SimCluster of its own? The
    primary is looked up among SimTracks only because, when it has no SimCluster,
    that is the one place it still appears.
    """
    section("4. Interior SimClusters")
    for key in KEYS:
        tree = trees[key]
        track_z = np.abs(np.concatenate(tree["SimTrack_boundaryPos_z"].array(library="np"))[
            np.concatenate(tree["SimTrack_crossedBoundary"].array(library="np")).astype(bool)])
        plane = front_face(track_z)
        print(f"\n{LABELS[key]}")

        def arr(branch):
            return tree[branch].array(library="np")

        sc_x, sc_y, sc_z = (arr(f"SimCluster_impactPoint_{a}") for a in "xyz")
        sc_cp = arr("SimCluster_CaloPartIdx")
        cp_energy, cp_pdg = arr("CaloPart_energy"), arr("CaloPart_pdgId")
        st_cross, st_pdg = arr("SimTrack_crossedBoundary"), arr("SimTrack_pdgId")
        st_z = arr("SimTrack_boundaryPos_z")
        st_pt, st_eta = arr("SimTrack_boundaryMomentum_pt"), arr("SimTrack_boundaryMomentum_eta")
        st_sc = arr("SimTrack_SimClusterIdx")

        n_cp = all_interior = parent_crossed = parent_kept = 0
        for ev in range(tree.num_entries):
            az = np.abs(sc_z[ev].astype(np.float64))
            faces, _ = classify_faces(az, radius(sc_x[ev], sc_y[ev]), key)
            crossed = st_cross[ev].astype(bool)
            for cp in range(len(cp_energy[ev])):
                n_cp += 1
                mine = sc_cp[ev] == cp
                if not mine.any() or not (faces[mine] == "interior").all():
                    continue
                all_interior += 1
                side = np.sign(sc_z[ev][mine][0])
                pmag = st_pt[ev] * np.cosh(st_eta[ev].astype(np.float64))
                parent = (crossed & (st_pdg[ev] == cp_pdg[ev][cp])
                          & front_face_mask(np.abs(st_z[ev]), plane)
                          & (np.sign(st_z[ev]) == side)
                          & (pmag > 0.5 * cp_energy[ev][cp]))
                if parent.any():
                    parent_crossed += 1
                    parent_kept += int((st_sc[ev][parent] >= 0).any())

        print(f"  CaloParticles whose every SimCluster is interior: {all_interior:,} "
              f"of {n_cp:,}")
        if all_interior:
            print(f"    primary crossed the front face (>50% of its energy): "
                  f"{parent_crossed:,}")
            print(f"    ...and that crossing became a SimCluster: {parent_kept:,}")


def figure_envelope(data, surfaces, material, tag, figdir, energy, dash="-"):
    """The envelope against the material it wraps, in the (|z|, r) plane.

    Two panels because the result lives at two scales: the 28.65 cm standoff is
    obvious in the full view and the 3-4 cm bore clearance is invisible there.

    Hand-rolled rather than routed through `hep_style.overlay` — this is a 2-D
    density with annotation, and `hep_style` deliberately offers only step
    histograms and a two-series ratio.
    """
    fig, (full, zoom) = panel_figure(2, figsize=(13.5, 5.6))

    for ax, (zlo, zhi), (rlo, rhi), title in (
        (full, (270, 525), (0, 170), "the boundary envelopes and the material"),
        # Wide enough to hold both bore steps and the radial jump between them:
        # the staircase is the shape claim, and one step alone does not show it.
        (zoom, (315, 412), (24, 40), "the bore staircase, and its clearance"),
    ):
        span = ((material["abs_z"] >= zlo) & (material["abs_z"] <= zhi)
                & (material["r"] >= rlo) & (material["r"] <= rhi))
        if span.any():
            ax.hist2d(material["abs_z"][span], material["r"][span], bins=(160, 120),
                      range=[[zlo, zhi], [rlo, rhi]], norm=mcolors.LogNorm(),
                      cmap="Greys", rasterized=True)

        for key in KEYS:
            d = data[key]
            plane = front_face(d["abs_z"])
            on_face = front_face_mask(d["abs_z"], plane)
            face_r = d["r"][on_face]
            # The front face, drawn over the radial extent it actually covers.
            ax.plot([plane, plane], [face_r.min(), face_r.max()],
                    color=COLORS[key], linewidth=2.2, solid_capstyle="butt", linestyle=dash,
                    label=series_label(key, "legacy", energy) if ax is full else None)
            for surface_r in surfaces[key]:
                on = ~on_face & (np.abs(d["r"] - surface_r) < SURFACE_TOL)
                if not on.any():
                    continue
                ax.plot([d["abs_z"][on].min(), d["abs_z"][on].max()],
                        [surface_r, surface_r], color=COLORS[key], linewidth=2.2,
                        solid_capstyle="butt", linestyle=dash)
            off = ~on_face
            ax.plot(d["abs_z"][off], d["r"][off], "o", color=COLORS[key],
                    markersize=2.6, linestyle="none", markeredgewidth=0, alpha=0.75)

        ax.set_xlim(zlo, zhi)
        ax.set_ylim(rlo, rhi)
        ax.set_xlabel("|z|  [cm]\n(distance along the beam axis)", linespacing=1.5)
        ax.set_ylabel("r  [cm]\n(distance from the beam axis)", linespacing=1.5)
        ax.set_title(title, pad=12)

    # The standoffs, drawn rather than left to the caption -- they are the result,
    # and both are measured here rather than written in.
    material_z = float(material["abs_z"].min())
    for key, height in (("fc0", 15.0), ("fc1", 150.0)):
        plane = front_face(data[key]["abs_z"])
        gap = material_z - plane
        # A span arrow needs room to read. FineCalo's 1.1 cm is a quarter of one
        # pixel on a 255 cm axis and renders as a scribble, so a gap that small
        # is labelled with a leader instead of drawn to scale.
        if gap > 0.04 * (525 - 270):
            full.annotate("", xy=(plane, height), xytext=(material_z, height),
                          arrowprops=dict(arrowstyle="<->", color=COLORS[key],
                                          linewidth=1.2, shrinkA=0, shrinkB=0))
            full.text(0.5 * (plane + material_z), height + 3.5,
                      f"{gap:.2f} cm", ha="center", fontsize=10, color=COLORS[key])
        else:
            full.annotate(f"flush, {gap:.2f} cm", xy=(plane, height),
                          xytext=(plane + 28, height + 12), fontsize=10,
                          color=COLORS[key],
                          arrowprops=dict(arrowstyle="->", color=COLORS[key],
                                          linewidth=1.0))

    # The bore clearance, on the panel where it is legible.
    d1 = data["fc1"]
    plane1 = front_face(d1["abs_z"])
    on_face1 = front_face_mask(d1["abs_z"], plane1)
    inner = [s for s in surfaces["fc1"] if s < np.median(d1["r"][on_face1])]
    if inner:
        bore = min(inner)
        on = ~on_face1 & (np.abs(d1["r"] - bore) < SURFACE_TOL)
        span = ((material["abs_z"] >= d1["abs_z"][on].min())
                & (material["abs_z"] <= d1["abs_z"][on].max()))
        edge = float(material["r"][span].min())
        mid = float(np.median(d1["abs_z"][on]))
        zoom.annotate("", xy=(mid, bore), xytext=(mid, edge),
                      arrowprops=dict(arrowstyle="<->", color=INK,
                                      linewidth=1.2, shrinkA=0, shrinkB=0))
        zoom.text(mid + 1.5, 0.5 * (bore + edge), f"{edge - bore:.2f} cm",
                  ha="left", va="center", fontsize=10, color=INK)

    full.legend(loc="upper right")
    full.annotate("grey: sensitive material (RecHits)", xy=(0.97, 0.80),
                  xycoords="axes fraction", ha="right", fontsize=10, color=INK_2)
    zoom.annotate("markers: crossings off the front face", xy=(0.5, 0.06),
                  xycoords="axes fraction", ha="center", fontsize=10, color=INK_2)
    return save(fig, figdir / f"boundary_envelope_{tag}.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="+", type=Path,
                        help="one nanoML file (its partner is found) or two")
    parser.add_argument("--tag", help="filename tag; default: the sample's")
    parser.add_argument("--material-events", type=int, default=MATERIAL_EVENTS,
                        help="events used for the RecHit material map")
    args = parser.parse_args()
    if len(args.files) > 2:
        parser.error("one or two files")

    sample = samples.load_sample(*args.files, tag=args.tag)
    trees = {key: sample.tree(key) for key in KEYS}
    figdir = style_dir(trees["fc0"], "envelope")
    for key in KEYS:
        print(f"{key}: {sample.files[key].name}")
    print(f"figures: {figdir}")

    data = {key: load(tree) for key, tree in trees.items()}
    section("0. The sample")
    for key in KEYS:
        print(f"{LABELS[key]:9s}: {len(data[key]['abs_z']):>8,} SimCluster crossings over "
              f"{trees[key].num_entries:,} events")

    tracks = {key: crossing_tracks(trees[key]) for key in KEYS}
    surfaces = calibrate(data, tracks)
    report_envelope(tracks, surfaces)

    material = material_map(trees["fc0"], args.material_events)
    report_standoff(data, tracks, surfaces, material)
    report_interior(trees)

    use_house_style()
    print("\nwrote", figure_envelope(data, surfaces, material, sample.tag, figdir, sample.energy,
                                     energy_linestyle(sample.energy)))


if __name__ == "__main__":
    main()
