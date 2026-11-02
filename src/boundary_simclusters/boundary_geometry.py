"""Where each config's boundary is, and which face of it a crossing landed on.

The boundary is **not a plane**. It is the closed surface of a 3-D region — a
thick annular disc with a conical hole punched down the middle for the beampipe,
the hole widening in coarse steps with depth. A boundary SimCluster is written
when a track crosses into that region, and the crossing can land on any face:

- the **front face**, a flat annulus at fixed |z| — 99.9% of crossings;
- the **inner bore**, a staircase of cylinders at fixed r going deeper;
- the **outer wall**, a cylinder at fixed large r.

The faces are corners of one shape, not separate features: FineCalo's front face
has inner edge r = 28.500 cm, exactly its first bore radius, and that bore begins
at |z| = 321.6 cm, just past the face at 321.047. Standard's front face has outer
edge r = 123.28 cm, matching its wall at 123.299, which runs up to its face at
293.5.

Two descriptions of the same geometry, and they agree. The *calibrated* surfaces
(`calibrate_surfaces`) are measured from the data: a constant-radius search, so
blind to every cone. The *source* envelopes (`SOURCE_ENVELOPES`, `ENVELOPE`) are
the two Geant4 volumes the mechanisms test -- `CALO` for Standard, `HGCal` for
FineCalo -- read from CMSSW's geometry XML. Every calibrated radius lands on a
source face to float32 precision, and the source adds what calibration cannot
reach: sloped faces, the risers between bore steps, and telling a point *inside*
a volume from one outside it.

A polycone's faces are all line segments in the (|z|, r) half-plane, so one
point-to-segment distance (`face_distances`) classifies planes, cylinders and
cones alike.

This module is imported, never run. It exists because the "is this crossing on
the plane" test was written out three times across the figure scripts with the
same magic 0.5, and the surface question needs one answer for the whole repo.

Nothing here reads a branch or opens a file — callers pass arrays, so the same
helpers serve every script.
"""

import numpy as np

# Half the 1 cm bin width of the boundary-z figures: wide enough to absorb
# float32 noise on a fixed surface, far narrower than the 27.5 cm between the two
# configs' front faces.
FRONT_FACE_TOL = 0.5

# A calibrated surface is exact to ~1e-3 cm (the spread within a group is
# 2-9e-4 cm across every local file), so this is ~50x the noise and still ~60x
# narrower than the smallest gap between two surfaces (28.50 -> 31.36).
SURFACE_TOL = 0.05

# A radius must be shared by at least this many crossings *and* this share of the
# non-front-face population to count as a surface.
MIN_SURFACE_COUNT = 3
MIN_SURFACE_SHARE = 0.01

# ...and every member must sit inside a window this wide. This is the guard that
# actually separates geometry from coincidence, and it is not optional: the
# photon fc0 sample has 2072 non-front-face crossings — the conversion
# vertices, a genuinely *continuous* distribution in r — and counting alone
# promoted five slices of it to surfaces. A real surface is an exactly repeated
# radius: the widest calibrated group in any local file spans 0.0023 cm, while a
# continuous population cannot put 1% of itself inside a 0.01 cm window.
MAX_SURFACE_SPREAD = 0.01

# Linkage distance for grouping radii before the spread test. Well below
# MAX_SURFACE_SPREAD so a real surface forms one group, and far below the
# spacing of anything else.
SURFACE_LINKAGE = 0.002

# ...and the group must hold this share of everything within SURFACE_TOL of its
# centre. The count and spread tests alone are fooled by a *fuzzy* population:
# single-linkage carves a continuum into slices, and on a small sample one slice
# can clear a floor of three or four. fc0 has such a population at r ~ 8.007 cm
# (spread 2e-2 cm, 100x a real surface) which was promoted to a surface on the
# 100-event pair and nowhere else. A real surface holds essentially everything
# near it -- measured concentration is 1.00 for all six real surfaces and
# 0.20-0.36 for the artefacts.
MIN_SURFACE_CONCENTRATION = 0.5

# The two Geant4 volumes `CaloTrkProcessing` tests, as (|z|, rMin, rMax) polycone
# planes in cm, positive-z half, exactly as the XML lists them. Read from CMSSW
# source by a separate session, not from this repo; the data then confirmed them.
#
# fc1 is the `HGCal` polycone, Geometry/HGCalCommonData/data/hgcal/v17/hgcal.xml
# with caloBase/2030/v7 constants. rMax is interpolated linearly between planes,
# so the outer wall is a chain of cones, not a staircase.
#
# fc0 is the `CALO` polycone, Geometry/CMSCommonData/data/caloBase.xml. It is the
# whole calorimeter system minus the tracker cavity, symmetric about z = 0, so it
# has no face at its first plane. Its forward bore is two cones of ~8 cm radius,
# ending at 526.5 cm in a step out to 272.5 cm (corrected 2026-11-01 from 554.1 cm).
#
# The version is pinned by the data, not assumed: all four calibrated fc1 bore
# radii and every bore-step |z| agree with v7 to float32 precision, and caloBase
# v8 would put the first step at 369.58 rather than the measured 366.51. An
# earlier copy of this geometry (31.4 / 44.8 / 154.0 / 438.8 / 262.0) mixed in
# v8 and v19 values; its cones missed the crossings by ~1.3 cm and needed a
# 2.5 cm tolerance to hide it.
#
# `names` label the faces `envelope_faces` derives, in its order: front cap,
# inner chain by z, rear cap, outer chain by z.
SOURCE_ENVELOPES = {
    "fc0": {
        "volume": "CALO",
        "mirrored": True,
        "planes": [(0.0, 123.3, 295.0), (293.5, 123.3, 295.0),
                   (293.5, 7.40, 295.0), (318.0, 8.00, 295.0),
                   (526.5, 8.93, 295.0), (526.5, 272.5, 295.0),
                   (554.1, 272.5, 295.0)],
        "names": ["tracker cavity wall", "tracker-cavity step",
                  "forward bore 1", "forward bore 2", "rear step", "rear bore",
                  "forward cap", "outer wall"],
    },
    "fc1": {
        "volume": "HGCal",
        "mirrored": False,
        "planes": [(321.05, 28.50, 152.33), (366.51, 28.50, 168.00),
                   (366.51, 31.36, 168.00), (387.75, 31.36, 175.32),
                   (410.73, 31.36, 204.43), (410.73, 38.36, 204.43),
                   (437.87, 38.36, 238.81), (437.87, 44.84, 238.81),
                   (456.20, 44.84, 263.45), (505.83, 44.84, 263.45),
                   (505.83, 44.84, 248.50), (523.70, 44.84, 248.50)],
        "names": ["front face",
                  "inner bore 1", "bore step 1", "inner bore 2", "bore step 2",
                  "inner bore 3", "bore step 3", "inner bore 4",
                  "rear face",
                  "outer cone 1", "outer cone 2", "outer cone 3",
                  "outer wall", "outer step", "outer rear wall"],
    },
}

# The face each volume's forward-going tracks overwhelmingly cross. Different
# surfaces of different volumes, which is the point: CALO's is the tracker-cavity
# step, 28.65 cm upstream of any calorimeter material; HGCal's is the calorimeter
# front face itself.
MAIN_FACE = {"fc0": "tracker-cavity step", "fc1": "front face"}

# Distance within which a crossing counts as *on* a source face. The same 0.05 cm
# as a calibrated surface, and it needs no special case for cones: crossings sit
# on the v7 outer cone to ~5e-3 cm and on the CALO forward bore to ~1e-3.
SOURCE_FACE_TOL = SURFACE_TOL

# Two consecutive segments of a chain whose slopes differ by less than this are
# one face. The XML repeats planes at every rMin/rMax change, and rounding to
# 0.01 cm moves a slope by ~1e-4; the smallest real slope change (HGCal's outer
# cone, 1.2668 -> 1.3443) is far above it.
COLLINEAR_TOL = 1e-3

C_CM_PER_NS = 29.9792458


def radius(x, y):
    """Cylindrical radius of a crossing, in float64.

    The surfaces are sharp to ~1e-3 cm and the branches are float32; promoting
    before the hypot keeps the calibration from smearing a surface into a
    distribution.
    """
    return np.hypot(np.asarray(x, dtype=np.float64),
                    np.asarray(y, dtype=np.float64))


def front_face(abs_z):
    """The |z| of this config's front face: the median of every crossing.

    Read from the data rather than hard-coded, because it is the thing under
    measurement — Standard lands at 293.5 cm and FineCalo at 321.047 cm, and a
    new sample must be allowed to say something else.
    """
    return float(np.median(abs_z))


def front_face_mask(abs_z, plane=None):
    """True where a crossing sits on the front face."""
    plane = front_face(abs_z) if plane is None else plane
    return np.abs(np.asarray(abs_z) - plane) < FRONT_FACE_TOL


def beta(x, y, z, t_seconds):
    """Path length from the origin over c*t — 1 for a prompt particle.

    These are gun events, so the origin is the right reference: a particle that
    flew straight from the interaction point at ~c arrives with beta ~ 1, while
    one that was slow (a sub-MeV neutron) or took a scattered path arrives late
    and reads low. It is a promptness proxy, not a velocity: a fast particle on a
    crooked path also reads below 1.
    """
    t_ns = np.asarray(t_seconds, dtype=np.float64) * 1e9
    distance = np.sqrt(np.asarray(x, dtype=np.float64) ** 2
                       + np.asarray(y, dtype=np.float64) ** 2
                       + np.asarray(z, dtype=np.float64) ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(t_ns > 0, distance / (C_CM_PER_NS * t_ns), np.nan)


def kinetic_energy(energy, pmag):
    """`boundaryEnergy` minus the rest mass it includes — see docs/data_handling.md.

    E - sqrt(E^2 - P^2). The clip is float32 hygiene: for a photon E and P are
    the same number to within rounding and the difference can go slightly
    negative.
    """
    energy = np.asarray(energy, dtype=np.float64)
    pmag = np.asarray(pmag, dtype=np.float64)
    return energy - np.sqrt(np.clip(energy ** 2 - pmag ** 2, 0.0, None))


def calibrate_surfaces(r_off, min_count=MIN_SURFACE_COUNT,
                       min_share=MIN_SURFACE_SHARE):
    """Radii shared by enough non-front-face crossings to be a real surface.

    `r_off` is the radius of every crossing that is *not* on the front face.
    Returns the surface radii, ascending.

    A real surface shows up as a radius repeated to within float32 noise by many
    crossings while z varies freely along it. A continuous population — photon
    conversion vertices, say — also puts many crossings at similar radii, so
    counting is not enough: the group must also be *narrow*. Both the floor and
    the spread test have to be cleared.
    """
    r_off = np.asarray(r_off, dtype=np.float64)
    if r_off.size == 0:
        return []
    floor = max(min_count, int(np.ceil(min_share * r_off.size)))

    # Single-linkage grouping on the sorted radii: consecutive values further
    # apart than SURFACE_LINKAGE start a new group. Rounding to a fixed grid
    # instead would split a surface that straddles a bin edge.
    ordered = np.sort(r_off)
    breaks = np.flatnonzero(np.diff(ordered) > SURFACE_LINKAGE) + 1
    surfaces = []
    for group in np.split(ordered, breaks):
        if group.size < floor or group[-1] - group[0] > MAX_SURFACE_SPREAD:
            continue
        # Is this a surface, or one slice of something continuous? A surface
        # holds nearly everything within SURFACE_TOL of it; a slice does not.
        centre = float(np.median(group))
        nearby = int(np.count_nonzero(np.abs(r_off - centre) < SURFACE_TOL))
        if group.size / max(nearby, 1) < MIN_SURFACE_CONCENTRATION:
            continue
        surfaces.append(centre)
    return sorted(surfaces)


def surface_name(surface_r, front_face_r):
    """`inner_bore` or `outer_wall`, by which side of the front face it is on.

    `front_face_r` is the median radius of the front-face crossings: a surface
    inside that is part of the beampipe hole, one outside is the envelope's
    outer skin.
    """
    kind = "inner_bore" if surface_r < front_face_r else "outer_wall"
    return f"{kind}@{surface_r:.2f}"


def _chain(points):
    """A polyline with repeated points dropped and collinear segments merged."""
    kept = [points[0]]
    for point in points[1:]:
        if point == kept[-1]:
            continue
        if len(kept) >= 2:
            (z0, r0), (z1, r1) = kept[-2], kept[-1]
            z2, r2 = point
            if z1 > z0 and z2 > z1 and abs((r1 - r0) / (z1 - z0)
                                           - (r2 - r1) / (z2 - z1)) < COLLINEAR_TOL:
                kept[-1] = point
                continue
        kept.append(point)
    return list(zip(kept, kept[1:]))


def _unit(nz, nr):
    norm = np.hypot(nz, nr)
    return nz / norm, nr / norm


def envelope_faces(config):
    """Every face of a source envelope as (name, z0, r0, z1, r1, nz, nr).

    The segment in the (|z|, r) half-plane plus its **outward** unit normal. The
    inner chain (rMin along z) bounds the volume from below in r, so its outward
    normal points to smaller r -- or along z for a step; the outer chain (rMax)
    is the mirror case. Deriving the normal from which chain a face belongs to is
    what gets CALO right: that volume is a shell whose centroid lies outside it,
    so a centroid heuristic flips every normal.
    """
    envelope = SOURCE_ENVELOPES[config]
    planes = envelope["planes"]
    first, last = planes[0], planes[-1]
    segments = []
    if not envelope["mirrored"]:
        segments.append(((first[0], first[1]), (first[0], first[2]), (-1.0, 0.0)))
    for a, b in _chain([(z, rmin) for z, rmin, _ in planes]):
        segments.append((a, b, _unit(b[1] - a[1], -(b[0] - a[0]))))
    segments.append(((last[0], last[1]), (last[0], last[2]), (1.0, 0.0)))
    for a, b in _chain([(z, rmax) for z, _, rmax in planes]):
        segments.append((a, b, _unit(-(b[1] - a[1]), b[0] - a[0])))

    names = envelope["names"]
    if len(names) != len(segments):
        raise ValueError(f"{config}: {len(segments)} faces but {len(names)} names")
    return [(name, a[0], a[1], b[0], b[1], nz, nr)
            for name, (a, b, (nz, nr)) in zip(names, segments)]


ENVELOPE = {key: envelope_faces(key) for key in SOURCE_ENVELOPES}

# A face that spans a range of |z| is a side wall -- a bore, a cone, the outer
# wall -- and a crossing on it entered or left through the *side* of the volume.
# A face at constant |z| is a transverse disc: the front and rear faces, the bore
# steps, the outer step. This is what "sideways" means in the census tables.
#
# Deliberately not "the outward normal is mostly radial": HGCal's outer cones 2
# and 3 have slope ~1.27, so that test files two thirds of the outer wall with
# the discs, and the front face of CALO's forward bore with the side walls.
SIDE_FACES = {key: frozenset(name for name, z0, _r0, z1, *_ in faces if z0 != z1)
              for key, faces in ENVELOPE.items()}


def face_distances(abs_z, r, faces):
    """Perpendicular distance from each crossing to each face segment, in cm."""
    p = np.stack([np.asarray(abs_z, float), np.asarray(r, float)], axis=1)
    out = np.empty((len(p), len(faces)))
    for i, (_, z0, r0, z1, r1, _nz, _nr) in enumerate(faces):
        a = np.array([z0, r0])
        ab = np.array([z1, r1]) - a
        t = np.clip((p - a) @ ab / (ab @ ab), 0.0, 1.0)
        out[:, i] = np.linalg.norm(p - (a + t[:, None] * ab), axis=1)
    return out


def classify_faces(abs_z, r, config, tol=SOURCE_FACE_TOL):
    """Nearest source face per crossing, with its distance.

    Beyond `tol` of every face a crossing is `interior` if it lies inside the
    volume -- no surface was crossed there at all -- and
    `exterior` otherwise.
    """
    faces = ENVELOPE[config]
    d = face_distances(abs_z, r, faces)
    best = d.argmin(axis=1)
    dist = d[np.arange(len(best)), best]
    names = np.array([f[0] for f in faces], dtype=object)[best]
    far = dist > tol
    names[far] = np.where(source_inside(np.asarray(abs_z)[far], np.asarray(r)[far],
                                        config), "interior", "exterior")
    return names, dist


def outward_normal(names, config):
    """The outward unit normal (nz, nr) of the face each crossing sits on.

    Zero for `interior` and `exterior`, which sit on no face.
    """
    nz = np.zeros(len(names))
    nr = np.zeros(len(names))
    for name, *_, a, b in ENVELOPE[config]:
        m = names == name
        nz[m], nr[m] = a, b
    return nz, nr


def source_inside(abs_z, r, config):
    """True where a point lies within a source envelope's volume.

    Radial limits are interpolated along each non-degenerate plane pair, which is
    how Geant4 builds a polycone. Points exactly on a riser are ambiguous, but
    those are on a face and never reach this test in `classify_faces`.
    """
    abs_z, r = np.asarray(abs_z, dtype=np.float64), np.asarray(r, dtype=np.float64)
    planes = SOURCE_ENVELOPES[config]["planes"]
    inside = np.zeros(abs_z.shape, dtype=bool)
    for (za, rmin_a, rmax_a), (zb, rmin_b, rmax_b) in zip(planes, planes[1:]):
        if za == zb:
            continue
        in_span = (abs_z >= za) & (abs_z <= zb)
        frac = (abs_z - za) / (zb - za)
        rmin = rmin_a + frac * (rmin_b - rmin_a)
        rmax = rmax_a + frac * (rmax_b - rmax_a)
        inside |= in_span & (r > rmin) & (r < rmax)
    return inside


# --- crossings, labelled by face -------------------------------------------------


def face_groups(config):
    """{face name: group} for one config's volume, from each face's geometry."""
    groups = {}
    for name, z0, _r0, z1, _r1, _nz, nr in ENVELOPE[config]:
        if name == MAIN_FACE[config]:
            groups[name] = "main face"
        elif z0 != z1:
            groups[name] = "inner side" if nr < 0 else "outer side"
        else:
            groups[name] = "rear and steps"
    return groups


def label_crossings(d, config):
    """Face, group, class and entering flag for crossings `d` (a DataFrame with
    x, y, z, abs_z, r and the momentum's pt, eta, phi) of one config."""
    face, _ = classify_faces(d.abs_z.to_numpy(), d.r.to_numpy(), config)
    groups = face_groups(config)
    group = np.array([groups.get(f, f) for f in face], dtype=object)
    cls = np.where(face == MAIN_FACE[config], "main face",
                   np.where(np.isin(face, ("interior", "exterior")), face, "other faces"))
    # Momentum in the (|z|, r) half-plane; a negative component along the face's
    # outward normal is a crossing into the volume.
    pt, eta, phi = (d[c].to_numpy() for c in ("pt", "eta", "phi"))
    x, y, z, r = (d[c].to_numpy() for c in ("x", "y", "z", "r"))
    pz = pt * np.sinh(eta) * np.sign(z)
    pr = (pt * np.cos(phi) * x + pt * np.sin(phi) * y) / np.maximum(r, 1e-9)
    nz, nr = outward_normal(face, config)
    return d.assign(face=face, group=group, **{"class": cls},
                    entering=(pz * nz + pr * nr) < 0)
