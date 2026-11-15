"""Paired per-CaloParticle comparison for one matched fc0/fc1 pair.

Everything in `pipeline.make_figures` is a
file-level aggregate: two independent distributions shown side by side. This
script is the *paired* view — each fc0 CaloParticle against its own
counterpart in fc1, through the event mapping written by
the match step (`pipeline.match_events`). It also gives the eta-region and momentum-threshold breakdowns,
so that a new sample gets every number in one run.

CaloParticles are paired through `calopart_matching_<tag>.parquet`, also
written by the match step: their stored order is not the same in both files
(docs/data_handling.md), so index `cp` in fc0 is not always index `cp` in fc1.

Either SimCluster table: the boundary table (`BoundarySimCluster_*`) when both
files carry it, else the legacy one (`SimCluster_*`); `--table`
chooses. Files are read through `samples.py`, so one file is enough -- its
partner is found -- and the output is `paired_calopart_<tag>.parquet` on the
legacy table, `paired_calopart_<tag>_bound.parquet` on the boundary table.

Read `docs/data_handling.md` before quoting any energy: `boundaryEnergy` is
kinematic and includes rest mass, `recEnergy` is an exclusive deposit,
`sumHitEnergy` counts a hit for a SimCluster and for its ancestors, so it must
never be summed as a budget.

The summary step of scripts/analyze.py (`--steps summary`), `run(sample, table)`
for the tests. The pair needs its mappings first (`pipeline.match_events`).
"""


import numpy as np
import pandas as pd

from boundary_simclusters import output_paths
from boundary_simclusters import samples  # noqa: E402

KEYS = samples.KEYS
LABELS = samples.LABELS

# The gun band, from GenPart_eta. Outside it only secondaries can live, which
# is what makes the split in `eta_regions` worth making.
BAND = (1.70, 2.70)

# Below this many Standard objects the ratio is counting noise, not a result.
MIN_COUNT = 25

PARTICLE_TYPES = (("photon", 22), ("e-", 11), ("e+", -11))
THRESHOLDS = (0.0, 1.0, 10.0, 30.0)


def load(sample, table):
    """{config: (SimCluster frame, CaloParticle frame)} from `samples.py`."""
    return {k: (samples.simcluster_frame(sample.tree(k), table, k),
                samples.caloparticle_frame(sample.tree(k))) for k in KEYS}


def build_paired(data, caloparts):
    """One row per matched CaloParticle, both configs side by side.

    `calopart` is the fc0 index; the fc1 side is read at its own index,
    `calopart_fc1`, from the CaloParticle matching."""
    paired = caloparts[["fc0_row", "fc1_row", "fc0_cp", "fc1_cp"]].astype(int).rename(
        columns={"fc0_cp": "calopart", "fc1_cp": "calopart_fc1"})
    cp0 = data["fc0"][1].set_index(["event", "cp"])
    own = cp0.loc[list(zip(paired.fc0_row, paired.calopart))]
    paired["pdgId"] = own.pdg.to_numpy().astype(int)
    paired["eta"] = own.eta.to_numpy()
    paired["gun_energy"] = own.energy.to_numpy()
    for key, (row, cp) in (("fc0", ("fc0_row", "calopart")), ("fc1", ("fc1_row", "calopart_fc1"))):
        sums = data[key][0].groupby(["event", "cp"]).agg(
            n=("energy", "size"), sumE=("energy", "sum"), recE=("rec", "sum"))
        got = sums.reindex(list(zip(paired[row], paired[cp])))
        paired[f"n_{key}"] = got.n.fillna(0).to_numpy().astype(int)
        paired[f"sumE_{key}"] = got.sumE.fillna(0.0).to_numpy()
        paired[f"recE_{key}"] = got.recE.fillna(0.0).to_numpy()
    paired = paired[["fc0_row", "fc1_row", "calopart", "calopart_fc1", "pdgId", "eta",
                     "gun_energy", "n_fc0", "sumE_fc0", "recE_fc0",
                     "n_fc1", "sumE_fc1", "recE_fc1"]]
    paired["abs_eta"] = paired.eta.abs()
    paired["delta_rel"] = (paired.sumE_fc1 - paired.sumE_fc0) / paired.sumE_fc0
    return paired


def section(title):
    print(f"\n{title}\n{'-' * len(title)}")


def report_sample(paired, data, mapping, n_events):
    section("0. The sample")
    pdgs = sorted(paired.pdgId.unique().tolist())
    per_event = paired.groupby(["fc0_row"]).size()
    print(f"events in file            : {n_events}")
    print(f"events paired             : {len(mapping)}")
    print(f"CaloParticles compared    : {len(paired)}")
    print(f"CaloParticle pdgId        : {pdgs}")
    print(f"CaloParticles per event   : {sorted(per_event.unique().tolist())}")
    print(f"gun energy (median)       : {paired.gun_energy.median():.1f} GeV")
    print(f"|eta| range               : {paired.abs_eta.min():.3f} - {paired.abs_eta.max():.3f}")

    # CaloParticle identity across configs. Weak by
    # construction for a fixed-multiplicity gun, but it is the check that rules
    # out an ancestry-level explanation for everything below.
    pdgs = {k: data[k][1].groupby("event").pdg.apply(lambda s: sorted(s.tolist()))
            for k in KEYS}
    same_count = all(pdgs["fc0"].get(r0, []) == pdgs["fc1"].get(r1, [])
                     for r0, r1 in zip(mapping.fc0_row.values, mapping.fc1_row.values))
    print(f"CaloParticle count+pdgId identical in both configs: {same_count}")


def report_multiplicity(paired):
    section("1. Boundary SimClusters per CaloParticle")
    table = pd.DataFrame(
        {
            LABELS[k]: {
                "total": int(paired[f"n_{k}"].sum()),
                "mean": round(float(paired[f"n_{k}"].mean()), 2),
                "median": float(paired[f"n_{k}"].median()),
                "p95": float(paired[f"n_{k}"].quantile(0.95)),
                "max": int(paired[f"n_{k}"].max()),
                "with zero": int((paired[f"n_{k}"] == 0).sum()),
            }
            for k in KEYS
        }
    ).T
    print(table.to_string())
    print(f"\nratio of means   (FineCalo/Standard) : "
          f"{paired.n_fc1.mean() / paired.n_fc0.mean():.2f}x")
    print(f"ratio of medians (FineCalo/Standard) : "
          f"{paired.n_fc1.median() / paired.n_fc0.median():.2f}x")
    print(f"per-CaloParticle ratio, median       : "
          f"{(paired.n_fc1 / paired.n_fc0.clip(lower=1)).median():.2f}x")
    print(f"FineCalo has more, per CaloParticle  : "
          f"{100 * (paired.n_fc1 > paired.n_fc0).mean():.1f}%")


def report_energy(paired):
    section("2. Sigma boundaryEnergy per CaloParticle")
    table = pd.DataFrame(
        {
            LABELS[k]: {
                "mean (GeV)": round(float(paired[f"sumE_{k}"].mean()), 3),
                "median (GeV)": round(float(paired[f"sumE_{k}"].median()), 3),
                "1st pct (GeV)": round(float(paired[f"sumE_{k}"].quantile(0.01)), 3),
                "min (GeV)": round(float(paired[f"sumE_{k}"].min()), 3),
                "energy per SimCluster (GeV)": round(
                    float(paired[f"sumE_{k}"].sum() / paired[f"n_{k}"].sum()), 4
                ),
            }
            for k in KEYS
        }
    ).T
    print(table.to_string())

    gun = float(paired.gun_energy.median())
    for k in KEYS:
        low = 100 * float((paired[f"sumE_{k}"] < 0.75 * gun).mean())
        print(f"{LABELS[k]:9s}: {low:.2f}% of CaloParticles below 75% of the gun energy")

    print(f"\nrelative difference (FineCalo - Standard) / Standard:")
    print(f"  median : {100 * paired.delta_rel.median():+.2f}%")
    print(f"  mean   : {100 * paired.delta_rel.mean():+.2f}%")
    print(f"  FineCalo exceeds Standard for {100 * (paired.delta_rel > 0).mean():.1f}% "
          "of CaloParticles")

    # Does the deficit scale with how finely FineCalo split the shower?
    by_split = paired.groupby(
        pd.qcut(paired.n_fc1, 5, duplicates="drop"), observed=True
    ).agg(
        CaloParticles=("delta_rel", "size"),
        mean_n_fc1=("n_fc1", "mean"),
        median_deficit_pct=("delta_rel", lambda s: 100 * s.median()),
    )
    by_eta = paired.groupby(
        pd.cut(paired.abs_eta, [1.7, 2.0, 2.2, 2.4, 2.7]), observed=True
    ).agg(
        CaloParticles=("delta_rel", "size"),
        median_deficit_pct=("delta_rel", lambda s: 100 * s.median()),
    )
    print("\nby number of FineCalo SimClusters (how finely the shower was split):")
    print(by_split.round(2).to_string())
    print("\nby |eta| of the CaloParticle:")
    print(by_eta.round(2).to_string())
    print(f"\ncorrelation(relative difference, FineCalo SimCluster count): "
          f"{paired.delta_rel.corr(paired.n_fc1):+.3f}")

    # recEnergy is forced to agree: both configs share the same RecHits, and each hit's
    # energy is split with fractions summing to one. It is the control that says a
    # boundaryEnergy deficit is bookkeeping, not lost energy.
    section("3. Sigma recEnergy per CaloParticle (control)")
    for k in KEYS:
        print(f"{LABELS[k]:9s}: mean {paired[f'recE_{k}'].mean():.3f} GeV")
    print(f"difference (FineCalo - Standard): "
          f"{(paired.recE_fc1 - paired.recE_fc0).mean():+.4f} GeV")


def flat(data, key, mapping, field):
    """`field` of every SimCluster in this config's paired events."""
    frame = data[key][0]
    rows = mapping.fc0_row.values if key == "fc0" else mapping.fc1_row.values
    return frame[field][frame.event.isin(rows)].to_numpy()


def report_eta_regions(data, mapping):
    section("4. Angular coverage: FineCalo/Standard count ratio by region")
    eta = {k: np.abs(flat(data, k, mapping, "eta")) for k in KEYS}
    regions = {
        f"in band, |eta| {BAND[0]:.2f}-{BAND[1]:.2f}": lambda a: (a >= BAND[0]) & (a <= BAND[1]),
        f"outside, |eta| > {BAND[1]:.2f}": lambda a: a > BAND[1],
        f"outside, |eta| < {BAND[0]:.2f}": lambda a: a < BAND[0],
    }
    rows = []
    for name, mask in regions.items():
        counts = {k: int(mask(eta[k]).sum()) for k in KEYS}
        rows.append({
            "region": name,
            "Standard": counts["fc0"],
            "FineCalo": counts["fc1"],
            "ratio": round(counts["fc1"] / counts["fc0"], 1) if counts["fc0"] else np.nan,
        })
    print(pd.DataFrame(rows).set_index("region").to_string())


def report_type_thresholds(data, mapping):
    section("5. Particle type, and where the excess lives in momentum")
    pdg = {k: flat(data, k, mapping, "pdg") for k in KEYS}
    pmag = {k: flat(data, k, mapping, "pmag") for k in KEYS}

    shares = pd.DataFrame(
        {
            LABELS[k]: {
                **{name: round(100 * float((pdg[k] == code).mean()), 2)
                   for name, code in PARTICLE_TYPES},
                "other": round(100 * float(
                    (~np.isin(pdg[k], [c for _, c in PARTICLE_TYPES])).mean()), 2),
            }
            for k in KEYS
        }
    ).T
    print("type shares of boundary SimClusters (%):")
    print(shares.to_string())

    rows = []
    for cut in THRESHOLDS:
        row = {"threshold (GeV)": "no cut" if cut == 0 else f"> {cut:g}"}
        for name, code in PARTICLE_TYPES:
            above = {k: int(((pdg[k] == code) & (pmag[k] > cut)).sum()) for k in KEYS}
            row[name] = (round(above["fc1"] / above["fc0"], 2)
                         if above["fc0"] >= MIN_COUNT else np.nan)
        rows.append(row)
    print("\nFineCalo/Standard count ratio as the momentum threshold is raised")
    print(f"(blank where Standard keeps fewer than {MIN_COUNT} objects):")
    print(pd.DataFrame(rows).set_index("threshold (GeV)").to_string())


def run(sample, table=None):
    """The paired per-CaloParticle comparison of `sample` on `table` (default:
    boundary where both files carry it): prints the report, writes
    paired_calopart_<tag>[_bound].parquet, and returns that table."""
    table = table or sample.table
    if table not in sample.tables:
        raise SystemExit(f"{table} table not in both files ({sample.tables})")
    mapping, caloparts = sample.mappings()

    for key in KEYS:
        print(f"{key}: {sample.files[key].name}")
    print(f"table: {table} ({samples.TABLES[table]}_*)")
    print(f"pairing: event_matching_{sample.tag}.parquet, calopart_matching_{sample.tag}.parquet")

    data = load(sample, table)
    paired = build_paired(data, caloparts)

    report_sample(paired, data, mapping, sample.tree("fc0").num_entries)
    report_multiplicity(paired)
    report_energy(paired)
    report_eta_regions(data, mapping)
    report_type_thresholds(data, mapping)

    suffix = "_bound" if table == "boundary" else ""
    out = output_paths.PROCESSED / f"paired_calopart_{sample.tag}{suffix}.parquet"
    paired.to_parquet(out, index=False)
    print(f"\nwrote {out}  ({len(paired)} rows)")
    return paired
