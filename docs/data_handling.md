# Reading the nanoML files correctly

The branches in these files can mislead: some names promise more than the branch holds, and
the two files of a pair do not line up row by row. This page lists what to do and why. The
code already follows it; read it before writing anything new against the files. Branch by
branch reference: `docs/schema.md`.

## Pairing the two files

**Events are not row-aligned.** Row `i` in the Standard file and row `i` in the FineCalo file are
usually different events: comparing them directly silently compares unrelated events.

**`run`, `luminosityBlock` and `event` do not identify events across files.** 

**Pair events on the generator record.** The sorted (`pdgId`, `pt`, `eta`, `phi`) of the
GenParticles is bit-identical across a true pair and unique per event in these samples, because
the guns draw η and φ per event. A second method, a fingerprint of RecHit count and total energy,
agrees on every event and remains as a cross-check, and the pairing stops if the two disagree.

**Pair CaloParticles on exact (`pdgId`, `eta`, `phi`), never by index.** Their stored order can differ
between the files, and `CaloPart_trackId` swaps with them. 

**Empty branches.** `CaloPart_GenPartIdx` is −1 and `CaloPart_simEnergy` is 0 for every CaloParticle in
every file. Neither carries information, so there is no direct link from a CaloParticle to its GenPart.

## The SimCluster energies

They answer different questions, and mixing them biases any energy-weighted result.

| branch | what it is | can you sum it? |
|---|---|---|
| `SimCluster_boundaryEnergy` | energy the particle carried at the crossing: kinematic, not a deposit | yes on the legacy collection |
| `SimCluster_recEnergy` | reconstructed hit energy assigned to this SimCluster alone | yes |


## Positions and directions

**`SimCluster_impactPoint_eta` is the momentum direction, not the position angle.** For where a crossing is,
compute η from `impactPoint_x/y/z`; for where it was going, use the branch. `impactPoint_phi` behaves the same way.
 |


Nothing in nanoML records a particle's fate after the boundary. With `trackId` unique within an event, so a second
crossing has nowhere to be stored, these files cannot test whether particles leave and re-enter. That needs the
Geant4 history in RECO.


## Other branches

- **`SimCluster_boundaryPmag`, `SimCluster_boundaryPt`:** the momentum at the crossing, one value per SimCluster.
- **Hit sharing:** hits are not an observable in this analysis, but the association is fractional. `RecHitHGC_SimClusterNumMatch`
  gives how many SimClusters share a hit, and `RecHitHGC_SimCluster_MatchIdx` and `_MatchQual` hold their indices and energy
  fractions, flattened so that `len(MatchIdx) == sum(NumMatch)`. The fractions sum to 1 for every matched hit.
- **Parent links:** nanoML has none, so ancestry comes from RECO (`scripts/tools/extract_ancestry.py`).
