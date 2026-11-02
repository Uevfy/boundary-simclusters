# Reading the nanoML files correctly

The branches in these files can mislead: some names promise more than the branch holds, and
the two files of a pair do not line up row by row. This page lists what to do and why. The
code already follows it; read it before writing anything new against the files. Branch by
branch reference: `docs/schema.md`.

Numbers come from the files named in the README. Unless a section says otherwise, they are from
the 200 GeV electron pair.

## Pairing the two files

**Events are not row-aligned.** Row `i` in the Standard file and row `i` in the FineCalo file are
usually different events: comparing them directly silently compares unrelated events.

**`run`, `luminosityBlock` and `event` do not identify events across files.** Joining on `event`
lands on a different physical event in 973 of 1000 cases. They are most likely per-job counters.

**Pair events on the generator record.** The sorted (`pdgId`, `pt`, `eta`, `phi`) of the
GenParticles is bit-identical across a true pair and unique per event in these samples, because
the guns draw η and φ per event. A second method, a fingerprint of RecHit count and total energy,
agrees on every event and remains as a cross-check, and the pairing stops if the two disagree.

| sample | events paired, of 1000 | events with CaloParticles in another order |
|---|---|---|
| e⁻ 200 GeV | 995 | 0 |
| e⁻ 20 GeV | 998 | 0 |
| γ 200 GeV | 971 | 1 |
| γ 20 GeV | 992 | 2 |

Both methods leave the same events unpaired, so those events are absent from the other file
rather than missed by the matcher.

The cause is most likely multithreaded production: a 100-event pair produced with one thread lines
up perfectly on row order, on `event` and on the fingerprint. One small sample does not prove this,
but the contrast with every multithreaded pair makes coincidence unlikely.

The GEN key depends on the gun. NanoML stores η and φ coarsely, so one particle's (η, φ) can repeat:
it is unique in only 994 of 1000 events. The key is unique here because the gun fires two independent
particles. A one-particle gun at around 1000 events per file, or a gun with every variable fixed, would
repeat keys. Repeated keys are left unmatched and reported, and the fingerprint is the fallback.

**Pair CaloParticles on exact (`pdgId`, `eta`, `phi`), never by index.** Their stored order can differ
between the files (table above), and `CaloPart_trackId` swaps with them. Pairing by index matters in
practice: in the swapped γ 200 GeV event it reads as two 1:1 decompositions, while paired correctly one
photon has 2 Standard SimClusters against 1 FineCalo, and the other 1 against 2.

**Empty branches.** `CaloPart_GenPartIdx` is −1 and `CaloPart_simEnergy` is 0 for every CaloParticle in
every file. Neither carries information, so there is no direct link from a CaloParticle to its GenPart.

## The SimCluster energies

They answer different questions, and mixing them biases any energy-weighted result.

| branch | what it is | can you sum it? |
|---|---|---|
| `SimCluster_boundaryEnergy` | energy the particle carried at the crossing: kinematic, not a deposit | yes on the legacy collection |
| `SimCluster_recEnergy` | reconstructed hit energy assigned to this SimCluster alone | yes |


## Rest mass and `SimCluster_mass`

**`boundaryEnergy` includes the rest mass, so soft hadrons look hard.** Among the crossings off the front face:

| species | `boundaryEnergy` median | kinetic energy median |
|---|---|---|
| photon | 0.026 GeV | 0.026 GeV |
| neutron | 0.951 GeV | **0.011 GeV** |
| proton | 1.198 GeV | 0.260 GeV |

Any threshold, sum or ranking across species is partly a ranking of rest masses. To ask how energetic an object
was, use `boundary_geometry.kinetic_energy(E, P)`, which is `E − sqrt(E² − P²)`.

**`SimCluster_mass` is not a rest mass.** 23.7% of its values are negative, and it reads ±0.0625 and ±0.0884 GeV
for electrons. Its meaning is undetermined, so do not use it. Compute `sqrt(E² − p²)` or look the mass up from
`SimCluster_pdgId`.

## Positions and directions

**`SimCluster_impactPoint_eta` is the momentum direction, not the position angle.** For where a crossing is,
compute η from `impactPoint_x/y/z`; for where it was going, use the branch. `impactPoint_phi` behaves the same way.

| crossings | η from `x/y/z`, median | branch, median | median difference |
|---|---|---|---|
| on the front face | 2.265 | 2.216 | 0.014 |
| off the front face | 3.156 | 1.785 | **1.360** |

The η figures and η-band results use the branch. On the front face, which holds 99.9% of crossings in the electron
samples, the difference is very unlikely to move a result, but this was not verified sample by sample. Do not quote
those numbers as positions.

**`SimTrack_lastPos_*` is not where a track ended.** For a track that crossed the boundary it is a copy of
`boundaryPos`, identical to within 10⁻³ cm for all 361,550 crossing tracks in the FineCalo file. `SimCluster_lastPos_*`
behaves the same way.

| track | `boundaryPos` | `lastPos` |
|---|---|---|
| crossed the boundary | the crossing point | a copy of the crossing point |
| did not cross | 0 (unset) | a real endpoint |

Nothing in nanoML records a particle's fate after the boundary. With `trackId` unique within an event, so a second
crossing has nowhere to be stored, these files cannot test whether particles leave and re-enter. That needs the
Geant4 history in RECO.

## The two SimTrack momenta

| branch | what it is |
|---|---|
| `SimTrack_pt`, `_eta`, `_phi` | momentum **at production** |
| `SimTrack_boundaryMomentum_*` | momentum **when the boundary record was taken**, 0 if there is none |

Evidence: a gun electron has p = 200.00 GeV in a 200 GeV sample in both configs, while its `boundaryMomentum` is far
lower and differs by config. A photon that radiates nothing has the two equal to six digits. The legs of a conversion
have production momenta averaging 200.0001 GeV in total over 1035 CaloParticles.

- **Do not sum production momenta over a decomposition.** SimCluster tracks nest, so the sum double-counts ancestry: 321 GeV
  per 200 GeV electron CaloParticle. It is additive only for a conversion pair.
- **`boundaryMomentum` is a momentum and `boundaryEnergy` an energy including rest mass.** They differ by 0.511 MeV for an
  electron and 940 MeV for a neutron.

## Other branches

- **`SimCluster_boundaryPmag`, `SimCluster_boundaryPt`:** the momentum at the crossing, one value per SimCluster.
- **`SimCluster_boundaryP4`:** a single small float per SimCluster, not a four-vector. Its meaning is unclear, so do not rely on it.
- **Hit sharing:** hits are not an observable in this analysis, but the association is fractional. `RecHitHGC_SimClusterNumMatch`
  gives how many SimClusters share a hit, and `RecHitHGC_SimCluster_MatchIdx` and `_MatchQual` hold their indices and energy
  fractions, flattened so that `len(MatchIdx) == sum(NumMatch)`. The fractions sum to 1 for every matched hit.
- **Parent links:** nanoML has none, so ancestry comes from RECO (`scripts/tools/extract_ancestry.py`).
