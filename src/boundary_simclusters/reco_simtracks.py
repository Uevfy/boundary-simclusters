"""Read Geant4's SimTracks and SimVertices from a RECO file, without CMSSW.

The nano files have no parent-track branches (production vertices only for the
tracks with a TrackingParticle, `TrackingPart_Vtx_*`), and their
SimTrack table is a filtered subset (27823 of 31001 tracks in the γ 20 GeV
Standard sample). The RECO file keeps the full `g4SimHits` collections, and with
them the ancestry: a track's `ivert` indexes the vertex list, and that vertex's
`parentTrackId` is the track's parent (-1 for a primary).

uproot cannot read these collections itself: they are `std::vector<SimTrack>` /
`std::vector<SimVertex>` written **memberwise**. For each data member in
streamer order (base class `CoreSimTrack` / `CoreSimVertex` first), the values of
all N objects follow one another:

- basic members (int, unsigned, uchar) are packed big-endian arrays;
- object members (event id, Lorentz and 3-vectors) are written one object at a
  time, each with a 4-byte byte count (flag 0x40000000) and a 2-byte version,
  plus a 4-byte checksum when the version is 0, and possibly nested. The plain
  numbers always sit at the end of each object's payload, so we take its last
  k bytes and skip by the byte count.

The layouts below are the streamer infos of `io_v1::SimTrack` v3 and
`io_v1::SimVertex` v3 (CMSSW_20_0_0_pre1). Every parse asserts that it consumed
the event's bytes exactly, so a layout change fails loudly rather than
returning shifted values.

`crossedBoundary` is bit 2 of `trackInfo` (`CROSSED_BOUNDARY`): it reproduces
the nano `SimTrack_crossedBoundary` flag on all 27823 tracks of the sample
above.

This module is imported, never run.
"""

import struct

import numpy as np
import uproot

TRACKS = "SimTracks_g4SimHits__HLT."
VERTICES = "SimVertexs_g4SimHits__HLT."
CROSSED_BOUNDARY = 4


def _entries(tree, collection):
    """Raw bytes of each event's entry of one collection's branch."""
    branch = tree[collection + "obj"]
    for i in range(branch.num_baskets):
        basket = branch.basket(i)
        data, offsets = basket.data, basket.byte_offsets
        for start, stop in zip(offsets[:-1], offsets[1:]):
            yield bytes(data[start:stop])


def _header(raw):
    """(number of objects, offset of the first member) of a memberwise vector."""
    _, version = struct.unpack(">IH", raw[:6])
    if not version & 0x4000:
        raise ValueError("vector is not written memberwise")
    p = 6
    element_version = struct.unpack(">H", raw[p:p + 2])[0]
    p += 2
    if element_version == 0:
        p += 4                                  # checksum instead of a version
    n = struct.unpack(">i", raw[p:p + 4])[0]
    return n, p + 4


def _basic(raw, p, n, fmt):
    dtype = np.dtype(fmt)
    values = np.frombuffer(raw, dtype=dtype.newbyteorder(">"), count=n, offset=p)
    return values.astype(dtype), p + n * dtype.itemsize


def _objects(raw, p, n, fmt, k):
    """The last `k` values of type `fmt` of each of `n` objects."""
    dtype = np.dtype(fmt)
    size = dtype.itemsize * k
    out = np.empty((n, k), dtype=dtype)
    for i in range(n):
        end = p + 4 + (struct.unpack(">I", raw[p:p + 4])[0] & 0x3FFFFFFF)
        out[i] = np.frombuffer(raw[end - size:end], dtype=dtype.newbyteorder(">"))
        p = end
    return out, p


def parse_tracks(raw):
    """One event's SimTracks as a dict of arrays, one row per track."""
    n, p = _header(raw)
    t = {}
    t["eventId"], p = _objects(raw, p, n, "u4", 1)
    t["trackId"], p = _basic(raw, p, n, "u4")
    t["pdgId"], p = _basic(raw, p, n, "i4")
    t["p4"], p = _objects(raw, p, n, "f8", 4)               # px, py, pz, E at the vertex
    t["ivert"], p = _basic(raw, p, n, "i4")
    t["igenpart"], p = _basic(raw, p, n, "i4")
    t["tkposition"], p = _objects(raw, p, n, "f8", 3)
    t["tkmomentum"], p = _objects(raw, p, n, "f8", 4)
    t["idAtBoundary"], p = _basic(raw, p, n, "i4")
    t["posAtBoundary"], p = _objects(raw, p, n, "f4", 4)    # x, y, z, t
    t["momAtBoundary"], p = _objects(raw, p, n, "f4", 4)    # px, py, pz, E
    t["trackInfo"], p = _basic(raw, p, n, "u1")
    if p != len(raw):
        raise ValueError(f"SimTrack layout mismatch: parsed {p} of {len(raw)} bytes")
    return t


def parse_vertices(raw):
    """One event's SimVertices as a dict of arrays, one row per vertex."""
    n, p = _header(raw)
    v = {}
    v["eventId"], p = _objects(raw, p, n, "u4", 1)
    v["pos"], p = _objects(raw, p, n, "f8", 4)              # x, y, z, t
    v["parentTrackId"], p = _basic(raw, p, n, "i4")
    v["vtxId"], p = _basic(raw, p, n, "u4")
    v["procType"], p = _basic(raw, p, n, "u4")
    if p != len(raw):
        raise ValueError(f"SimVertex layout mismatch: parsed {p} of {len(raw)} bytes")
    return v


def read(path):
    """(tracks, vertices): one dict of arrays per event, in file order."""
    tree = uproot.open(path)["Events"]
    tracks = [parse_tracks(raw) for raw in _entries(tree, TRACKS)]
    vertices = [parse_vertices(raw) for raw in _entries(tree, VERTICES)]
    return tracks, vertices


def ancestry(tracks, vertices):
    """Per event, {trackId: parent trackId}; -1 for a primary."""
    return [dict(zip(t["trackId"].astype(int).tolist(),
                     v["parentTrackId"][t["ivert"]].astype(int).tolist()))
            for t, v in zip(tracks, vertices)]


def check_alignment(tracks, nano):
    """Raise unless every nano SimTrack sits in the same RECO event with the same
    pdgId: the nano file must be the nanoML of this RECO file, row for row.
    `nano` holds the arrays `SimTrack_trackId` and `SimTrack_pdgId`."""
    bad = []
    for i, t in enumerate(tracks):
        pdg = dict(zip(t["trackId"].tolist(), t["pdgId"].tolist()))
        if any(pdg.get(int(tid)) != int(p)
               for tid, p in zip(nano["SimTrack_trackId"][i], nano["SimTrack_pdgId"][i])):
            bad.append(i)
    if bad or len(tracks) != len(nano["SimTrack_trackId"]):
        raise SystemExit(f"RECO and nano rows do not line up ({len(bad)} events, "
                         f"first {bad[:5]}); is the nano made from this RECO file?")
