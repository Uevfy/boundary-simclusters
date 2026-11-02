# FineCalo vs. Standard CMSSW: boundary SimCluster comparison

A starting point for analysing nanoML ROOT files from HGCAL simulation. It reads the files,
pairs them across two configurations, compares their SimClusters and draws the standard
figures, and you can build new studies on top of it. The question it was built for is whether
two CMSSW mechanisms hand a later step the same inputs for identical simulated events: **do FineCalo (`useFineCalo=1`, "fc1") and the Standard mechanism
(`useFineCalo=0`, "fc0") give the same decomposition of a CaloParticle into boundary
SimClusters?**

The definition of a shower rooted at the calorimeter boundary follows K. Schaeuble et al.,
"Detector-aware target definitions for full-event particle reconstruction"
(arXiv:2606.30485). A CaloParticle groups everything from one primary particle, and can
hold several boundary SimClusters, for example a photon converting just before the
calorimeter.

Figures from this project are **not CMS-approved**.

## Where to start

1. Follow [Setup](#setup), get the ROOT files (see [Data](#data)), and run `pytest` to confirm
   your environment reproduces the stored numbers.
2. Read `docs/data_handling.md`. The branches have traps, such as events that do not line up
   between the files and an `eta` branch that is not a position.
3. Run `scripts/analyze.py` on one pair, then look through its report and figures.
4. Write your own study, as described in [Adding a study](#adding-a-study).

## Read this before comparing numbers

The two configs do not record on the same surface. Standard tests the Geant4 volume `CALO`
and FineCalo tests `HGCal`, whose front faces lie at |z| = 293.5 cm and 321.05 cm, with
material in between. Standard also flags a track whose first step is already inside `CALO`.
A raw fc0-vs-fc1 count or spectrum therefore mixes the mechanism with the surface. Energy is
not part of the classification: an electron loses a median of ~70% of its kinetic energy in
that gap.

## Setup

Python 3.10 or newer; the versions in `requirements.txt` are the tested ones. Nothing here
needs CMSSW or ROOT: `uproot` reads the files.

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -e .          # makes src/boundary_simclusters importable
pytest                    # skips every sample whose input files are absent
```

## Data

**The samples are not in this repository.** Without them the code can be read and installed, but
not run.
Simulation and production happen outside this repository, so nothing here can regenerate
them.

All samples are two-particle guns in |η| 1.70–2.70, 1000 events, seed 1, with
`CMSSW_20_0_0_pre1`, global tag `150X_mcRun4_realistic_v1`, HGCal v17 and caloBase 2030/v7. They come in pairs that differ only in the FineCalo
switches. The nanoML files go in the repo root, the ancestry files in `data/processed/`:

| tag | gun | nanoML (`bound_*`, both SimCluster collections) | RECO (only for ancestry) |
|---|---|---|---|
| `e200` | e⁻ 200 GeV | `bound_electron_fc{0,1}_200GeV_testNanoML.root` | `testRECO.root`, `FinetestRECO.root` |
| `e20` | e⁻ 20 GeV | `bound_electron_fc{0,1}_20GeV_testNanoML.root` | `22_fc{0,1}_20GeV_testRECO.root` |
| `g200` | γ 200 GeV | `bound_photon_fc{0,1}_200GeV_testNanoML.root` | `11_fc{0,1}_s1_testRECO.root` |
| `g20` | γ 20 GeV | `bound_photon_fc{0,1}_20GeV_testNanoML.root` | `photon_fc{0,1}_20GeV_testRECO.root` |

The golden test's fifth sample, `s1`, is `testNanoML.root` and `FinetestNanoML.root`: the
200 GeV electrons with only the legacy collection.

NanoML has no parent links, so the classification reads them from
`data/processed/ancestry_<tag>_<fc0|fc1>.parquet`, written once from RECO by
`scripts/tools/extract_ancestry.py`. Without those files a differing CaloParticle can only be
called `different tracks`, and the run warns.

Numeric filename prefixes are production labels, not pdgIds: `11_` is the photon gun and
`22_` the electron gun. The code reads particle, energy and config from the file contents,
never from the name.

## Running it

```bash
python scripts/analyze.py bound_photon_fc0_20GeV_testNanoML.root bound_photon_fc1_20GeV_testNanoML.root --tag g20
```

Give both files of a pair, in either order; with one, its partner is searched for in the same
folder. The default is the boundary collection (`BoundarySimCluster_*`), and `--table legacy`
reads `SimCluster_*`. `--steps` runs a subset and `--only` a subset of the figures. Everything
printed is saved as `report_<tag>[_bound].txt` beside the figures.

| `--steps` | module in `boundary_simclusters.pipeline` | what it does | writes |
|---|---|---|---|
| `match` | `match_events` | pairs events on the generator record, then CaloParticles on exact (pdgId, η, φ). The files are not row-aligned and their CaloParticle order can differ (`docs/data_handling.md`) | `data/processed/{event,calopart}_matching_<tag>.parquet` |
| `summary` | `paired_calopart_summary` | the paired per-CaloParticle comparison | `data/processed/paired_calopart_<tag>[_bound].parquet` |
| `figures` | `make_figures` | position, η, energy, pT, species, multiplicity and spectra | `outputs/<particle>/<energy>/<table>/` |
| `surfaces` | `classify_surfaces` | which face of its config's volume each crossing hit | figures, `data/processed/surface_<tag>[_bound].parquet` |
| `decompositions` | `classify_decompositions` | compares each CaloParticle's SimClusters by root track, one class each | `data/processed/decompositions_<tag>[_bound].parquet`, `decomposition_simclusters_<tag>[_bound].parquet` |
| `compare` | `make_figures`, `classify_surfaces` | the figures again with both collections overlaid | `outputs/<particle>/<energy>/comparisons/` |

### How SimClusters are compared

Across configs, within a matched CaloParticle, SimClusters are compared by the Geant4 track
that opened each. A track stored in both files of a pair carries the same ID, because both
configs simulate the same physics from the same seed; `id_agreement` repeats this check on
every run and stops if it fails. Track IDs restart per event, so they cannot pair events. Each
CaloParticle gets one class on the boundary collection:

| class | meaning |
|---|---|
| same tracks | both configs record the same root tracks |
| photon conversion | every replaced Standard root is a photon, and its FineCalo descendants come through its conversion |
| other replacement | a Standard root is replaced by its FineCalo descendants, not only through a conversion: an object lost in the gap |
| no descendant | a Standard root has no FineCalo descendant |
| mixed | a replacement and a Standard root with no descendant |
| unlinked | a FineCalo-only root has no Standard-only ancestor |

The legacy collection has its own classes (`LEGACY_CLASSES` in `classify_decompositions`).

## Studies and tools

Each takes nanoML files; the classes they read come from `analyze.py`.

| script | what it does | input |
|---|---|---|
| `studies/topology_breakdown.py` | the classes by species, energy and \|η\|, with binomial errors; figures vs \|η\| | `bound_*` files, after `analyze.py` |
| `studies/classify_boundary_surface.py` | each surface measured against the CMSSW source radii, its standoff from the sensitive material, the envelope figure | a pair |
| `studies/reco_crosscheck.py` | radiation lengths before and in the gap, from where gun photons convert | both RECO files and both nanoML files of a pair |
| `studies/check_boundary_nesting.py` | tests the boundary-SimCluster rule against full Geant4 ancestry | a RECO file and its nanoML |
| `studies/compare_collections.py` | boundary vs legacy collection on the same events | `bound_*` files, after `analyze.py` on both tables |
| `studies/event_scatter.py` | event display of one event in both configs; draws only what the branches say | a pair and `--row` |
| `studies/decomposition_example.py` | one CaloParticle's decomposition into SimClusters, per collection | a pair, after `analyze.py --steps match` |
| `studies/compare_productions.py` | one gun sample from two productions, event by event. It needs each production's ancestry for classes; the per-CaloParticle comparison uses counts only | two pairs, `--root <tree> --tag <tag>` |
| `tools/extract_ancestry.py` | writes the ancestry parquets from RECO | both nanoML and both RECO files of a pair |
| `tools/print_root_branches.py` | branch names and types, the source of `docs/schema.md` | nanoML files |

## Adding a study

Put a new script in `scripts/studies/` and keep it thin; anything other code will call belongs
in `src/boundary_simclusters/`. Read files through `samples` instead of opening branches by hand:

```python
from boundary_simclusters import samples

sample = samples.load_sample("bound_photon_fc0_20GeV_testNanoML.root",
                             "bound_photon_fc1_20GeV_testNanoML.root")
frame = samples.simcluster_frame(sample.tree("fc1"), "boundary", "fc1")   # one row per SimCluster
```

Compare across the two files only through the mappings from `analyze.py --steps match`
(`sample.mappings()`), never by row or index. For figures, use `hep_style` and save through
`output_paths`, which works out the folder from the data. After any change under `src/`, run
`pytest`. If a number moved on purpose, regenerate with `python tests/golden.py` and commit the new
`tests/golden.json` with the change (`tests/README.md`).

## Layout

```
src/boundary_simclusters/     the library, imported and never run
    pipeline/                 one module per step, each with run(sample, ...)
    samples.py                what a file is, its partner, the per-SimCluster readers
    matching.py               event and CaloParticle pairing
    boundary_geometry.py      both boundary volumes as CMSSW defines them, face labels
    reco_simtracks.py         Geant4 tracks and parent links from RECO, without CMSSW
    hep_style.py              the figure style
    output_paths.py           where each output belongs, read from the data
scripts/                      analyze.py, studies/, tools/
tests/                        the regression test (tests/README.md)
docs/                         see Documentation below
```

`outputs/` and `data/processed/` are created on the first run and stay out of git.

**Figures.** `hep_style.save` writes each figure as a PNG with its title and no experiment label.
A PDF version with the label "Private work (CMS simulation)" is commented out in `save`.

## Tests

`pytest` reruns the pipeline on the five samples and compares every number with
`tests/golden.json`, in about a minute. It proves the code still gives what it gave, not that
it is right. See `tests/README.md`.

## Documentation

- `docs/data_handling.md`: how to read the files correctly, including event pairing and what each
  energy branch means
- `docs/schema.md`: the nanoML branches, field by field
