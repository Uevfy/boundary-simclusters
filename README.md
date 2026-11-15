# FineCalo vs. Standard CMSSW: boundary SimCluster comparison

A starting point for analysing nanoML ROOT files from HGCAL simulation. It reads the files,
pairs them across two configurations, compares their SimClusters and draws the standard
figures, and you can build new studies on top of it. The question it was built for is whether
two CMSSW mechanisms hand a later step the same inputs for identical simulated events: **do FineCalo (`useFineCalo=1`, "fc1") and the Standard mechanism
(`useFineCalo=0`, "fc0") give the same decomposition of a CaloParticle into boundary
SimClusters?**

Figures from this project are **not CMS-approved**.

## Where to start

1. Follow [Setup](#setup), then produce a sample pair as described in [Data](#data).
2. Run `scripts/analyze.py` on one pair, then look through its report and figures.
3. Write your own study, as described in [Adding a study](#adding-a-study).

## Setup

Python 3.11 or newer (the pinned pandas and numpy need it); the versions in
`requirements.txt` are the tested ones. Nothing here needs CMSSW or ROOT: `uproot` reads
the files.

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -e .          # makes src/boundary_simclusters importable
```

## Data

**No samples come with this repository: you produce them yourself.** Simulation and
production run in CMSSW, outside this repository. Use these branches:

| repository | branch | what it holds |
|---|---|---|
| [Uevfy/production_tests](https://github.com/Uevfy/production_tests/tree/d121-production) | `d121-production` | `GSD_GUN.py`, `RECO.py`, `nanoML_cfg.py` with the BoundarySimCluster table, batch submission. Its README is the production recipe |
| [Uevfy/cmssw](https://github.com/Uevfy/cmssw/tree/thesis-production) | `thesis-production` | the pepr `CMSSW_20_0_0_pre1` used for the samples |
| [Uevfy/reco-prodtools](https://github.com/Uevfy/reco-prodtools/tree/thesis-production) | `thesis-production` | the fragments `GSD_GUN.py` and `RECO.py` import (GSD Standard and FineCalo, RECO) as used for the samples, and `produceSkeletons_D121.sh` |

Follow the setup in the `d121-production` README, but take these forks instead of the
`dgaytanv` ones it names:

```bash
git cms-merge-topic Uevfy:thesis-production
git clone -b thesis-production https://github.com/Uevfy/reco-prodtools.git reco_prodtools
git clone -b d121-production https://github.com/Uevfy/production_tests.git
```

A pair is the same gun run twice with the same seed, once per config (`useFineCalo=0` for
Standard, `useFineCalo=1` for FineCalo), each through GSD, RECO and nanoML. The original
results come from two-particle guns in |η| 1.70–2.70, 1000 events, seed 1: e⁻ and γ at 20 and
200 GeV, `CMSSW_20_0_0_pre1`, global tag `150X_mcRun4_realistic_v1`, HGCal v17 and caloBase
2030/v7. Those samples were made with `master-production`; `d121-production` gives the same
SimClusters event by event.

**Keep the RECO files.** NanoML has no parent links, so the classification reads them from
RECO, through `data/processed/ancestry_<tag>_<fc0|fc1>.parquet`:

```bash
python scripts/tools/extract_ancestry.py photon_20GeV_fc0_NanoML.root photon_20GeV_fc1_NanoML.root \
    photon_20GeV_fc0_RECO.root photon_20GeV_fc1_RECO.root
```

Without those files a CaloParticle whose root tracks differ can only be called
`different tracks`, and the run warns.

### Naming files and tags

The code reads nothing from a file name. Standard or FineCalo, particle and energy come from the
file contents, and each RECO file is tied to its nanoML file by matching rows, so any name ending
in `.root` works. A readable scheme, used in the examples here:

```
<particle>_<energy>_fc<0|1>_NanoML.root     photon_20GeV_fc0_NanoML.root
<particle>_<energy>_fc<0|1>_RECO.root       photon_20GeV_fc0_RECO.root
```

with `fc0` for Standard (`useFineCalo=0`) and `fc1` for FineCalo (`useFineCalo=1`). Two things to know:

- **The tag** is chosen for you. Every output of a pair carries a label, its tag, built from the
  contents as `<particle>_<energy>`: `ancestry_photon_20GeV_fc0.parquet`,
  `report_photon_20GeV_bound.txt`, `n_per_event_photon_20GeV_bound.png`. Every script picks the
  same one, so the steps find each other's outputs. Give `--tag` only to keep apart two samples
  of the same particle and energy, such as another seed or another production, and then give
  that same tag to every script of the pair.
- **The folder.** Given one file, `analyze.py` looks for its partner among the `.root` files in
  the same folder, by shared events. Two productions of the same gun and seed in one folder make
  that ambiguous, and it stops and asks for both files. Keep each production in its own folder.

## Running it

```bash
python scripts/analyze.py photon_20GeV_fc0_NanoML.root photon_20GeV_fc1_NanoML.root
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

sample = samples.load_sample("photon_20GeV_fc0_NanoML.root",
                             "photon_20GeV_fc1_NanoML.root")
frame = samples.simcluster_frame(sample.tree("fc1"), "boundary", "fc1")   # one row per SimCluster
```

Compare across the two files only through the mappings from `analyze.py --steps match`
(`sample.mappings()`), never by row or index. For figures, use `hep_style` and save through
`output_paths`, which works out the folder from the data. After any change under `src/`, rerun
`analyze.py` on a pair and compare its report with the one from before the change.

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
docs/                         see Documentation below
```

`outputs/` and `data/processed/` are created on the first run and stay out of git.

**Figures.** `hep_style.save` writes each figure as a PNG with its title and no experiment label.
A PDF version with the label "Private work (CMS simulation)" is commented out in `save`.

## Documentation

- `docs/data_handling.md`: how to read the files correctly, including event pairing and what each
  energy branch means
- `docs/schema.md`: the nanoML branches, field by field
