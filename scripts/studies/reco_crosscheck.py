"""The material before and between the two surfaces, from where each gun photon converts.

    python scripts/studies/reco_crosscheck.py photon_fc0_20GeV_testRECO.root \
        photon_fc1_20GeV_testRECO.root bound_photon_fc0_20GeV_testNanoML.root \
        bound_photon_fc1_20GeV_testNanoML.root

Takes both configs' RECO files and the nanoML made from each, in any order: RECO
and nanoML are told apart by their branches, the config is read from the nanoML,
and each RECO is paired with the nanoML whose rows it lines up with
(`reco_simtracks.check_alignment`). The probability to convert before Standard's
surface (|z| < 293.5 cm) measures the tracker's radiation lengths along that
|eta|; the probability to convert in the gap given the photon arrived unconverted
measures the gap's (P = 1 - exp(-7 X0 / 9)). Read from FineCalo's RECO, which
stores every track crossing 321.05 cm, so every conversion before it. (Vertices
of *stored* tracks are no material map: each config stores the tracks crossing
its own boundary.)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

from boundary_simclusters import samples
from boundary_simclusters.reco_simtracks import TRACKS, check_alignment, read

STD_FACE, FINE_FACE = 293.5, 321.05            # |z| of each config's main face [cm]
ETA_BINS = np.linspace(1.70, 2.70, 6)
CONVERSION = 14                                # Geant4 process subtype of SimVertex.procType


def sort_files(paths):
    """{config: (RECO path, nanoML path)} from four paths in any order."""
    reco, nano = [], {}
    for path in paths:
        keys = set(uproot.open(path)["Events"].keys(recursive=False))
        if TRACKS in keys:
            reco.append(path)
        elif "SimTrack_trackId" in keys:
            nano[samples.detect_config(samples.events(path))] = path
    if len(reco) != 2 or set(nano) != set(samples.KEYS):
        raise SystemExit("need two RECO files and the Standard and FineCalo nanoML made from them")
    out = {}
    for key, nano_path in nano.items():
        arrays = samples.events(nano_path).arrays(["SimTrack_trackId", "SimTrack_pdgId"], library="np")
        for path in reco:
            tracks, vertices = read(path)
            try:
                check_alignment(tracks, arrays)
            except SystemExit:
                continue
            out[key] = (path, nano_path, tracks, vertices)
            break
        else:
            raise SystemExit(f"no RECO file lines up with {nano_path.name}")
    return out


def gun_photon_conversions(tracks, vertices, nano_path):
    """One row per gun photon (a CaloParticle's track with pdgId 22): its |eta| and
    the |z| of its conversion vertex (inf if none is stored)."""
    calo = samples.events(nano_path)["CaloPart_trackId"].array(library="np")
    rows = []
    for event, (t, v) in enumerate(zip(tracks, vertices)):
        index = {k: i for i, k in enumerate(t["trackId"].astype(int))}
        for primary in calo[event].astype(int):
            i = index.get(primary)
            if i is None or t["pdgId"][i] != 22:
                continue
            px, py, pz = t["p4"][i, :3]
            converted = (v["parentTrackId"] == primary) & (v["procType"] == CONVERSION)
            z = np.abs(v["pos"][converted, 2]).min() if converted.any() else np.inf
            rows.append((abs(np.arcsinh(pz / np.hypot(px, py))), z))
    return pd.DataFrame(rows, columns=["eta", "z_conversion"])


def radiation_lengths(p):
    return -9 / 7 * np.log(1 - p)


def report_material(photons):
    print("\n1. Material before and between the surfaces, from where gun photons convert")
    print("---------------------------------------------------------------------------")
    if not len(photons):
        print("no gun photons in this sample")
        return
    rows = {}
    for interval, g in photons.groupby(pd.cut(photons.eta, ETA_BINS), observed=True):
        before = g.z_conversion < STD_FACE
        arrived = g[~before]
        in_gap = (arrived.z_conversion >= STD_FACE) & (arrived.z_conversion < FINE_FACE)
        pb, pg = before.mean(), in_gap.mean()
        rows[f"{interval.left:.1f}-{interval.right:.1f}"] = {
            "photons": len(g),
            "converted before 293.5 cm": f"{100 * pb:5.1f} +- {100 * np.sqrt(pb * (1 - pb) / len(g)):3.1f}%",
            "X0 upstream": round(radiation_lengths(pb), 2),
            "converted in the gap | arrived": f"{100 * pg:5.1f} +- {100 * np.sqrt(pg * (1 - pg) / len(arrived)):3.1f}%",
            "X0 in the gap": round(radiation_lengths(pg), 2)}
    print(pd.DataFrame(rows).T.to_string())
    overall = ((photons.z_conversion >= STD_FACE) & (photons.z_conversion < FINE_FACE)).mean()
    print(f"converted in the gap: {100 * overall:.1f}% of all gun photons "
          f"(compare: gap conversions + single descendants in the classification)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs=4, type=Path,
                    help="both RECO files and the nanoML made from each, any order")
    args = ap.parse_args()
    files = sort_files(args.files)
    for key in samples.KEYS:
        print(f"{samples.LABELS[key]:9s} {files[key][0].name} <-> {files[key][1].name}")
    report_material(gun_photon_conversions(*files["fc1"][2:], files["fc1"][1]))


if __name__ == "__main__":
    main()
