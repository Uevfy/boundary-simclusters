# ROOT/nanoML schema reference

All 414 branches of the `Events` tree of the `bound_*` nanoML files, grouped by
collection and by relation within each collection. Identical (name, type) in
both configs. The older files without the boundary collection (`testNanoML.root`,
`FinetestNanoML.root`, `11_*`, `22_*`, `photon_*`) carry the first 288 of them:
everything except the sections and rows marked *`bound_*` only*.
Regenerate/verify with `scripts/tools/print_root_branches.py`.

Types are as reported by `uproot`; `[]` means one jagged array entry per
object per event, no `[]` means one scalar value per event.

Collections **not** needed for the boundary-SimCluster question are grouped at
the bottom of this file, under "Reference only". Everything above that divider
is a collection this analysis actually reads; nothing is omitted, only
reordered.

---

## CaloPart (17) — CaloParticle. Paper concept: primary particle passed to GEANT4.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `CaloPart_pdgId` | `int32_t[]` | |
| Identity | `CaloPart_charge` | `int32_t[]` | |
| Identity | `CaloPart_GenPartIdx` | `int32_t[]` | index into `GenPart` — **−1 for every CaloParticle in every local file**, unusable (docs/data_handling.md) |
| Identity | `CaloPart_trackId` | `int32_t[]` | Geant4 track id: a track stored in both files of a pair has the same id in both (`classify_decompositions.id_agreement` checks this on every run); ids restart per event; its index can differ between the files |
| Identity | `CaloPart_eventId` | `int32_t[]` | |
| Identity | `CaloPart_rawEventId` | `int32_t[]` | |
| Identity | `CaloPart_bunchCrossing` | `int32_t[]` | |
| Identity | `CaloPart_nGenPart` | `int32_t[]` | |
| Identity | `CaloPart_nSimHit` | `int32_t[]` | |
| Identity | `CaloPart_nSimTrack` | `int32_t[]` | |
| Boundary | `CaloPart_crossedBoundary` | `bool[]` | |
| Kinematics | `CaloPart_energy` | `float[]` | |
| Kinematics | `CaloPart_eta` | `float[]` | |
| Kinematics | `CaloPart_phi` | `float[]` | |
| Kinematics | `CaloPart_pt` | `float[]` | |
| Kinematics | `CaloPart_mass` | `float[]` | |
| Kinematics | `CaloPart_simEnergy` | `float[]` | **0 for every CaloParticle in every local file** (docs/data_handling.md); `mass` unreliable if analogous to `SimCluster_mass` |

## GenPart (9) — GEN-level particle record, pre-simulation.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `GenPart_pdgId` | `int32_t[]` | |
| Identity | `GenPart_charge` | `int32_t[]` | |
| Identity | `GenPart_status` | `int32_t[]` | |
| Identity | `GenPart_statusFlags` | `uint16_t[]` | bitmask |
| Identity | `GenPart_genPartIdxMother` | `int16_t[]` | |
| Kinematics | `GenPart_eta` | `float[]` | |
| Kinematics | `GenPart_phi` | `float[]` | |
| Kinematics | `GenPart_pt` | `float[]` | |
| Kinematics | `GenPart_mass` | `float[]` | |

## GenVtx (4) — primary GEN vertex, one per event (not an array).

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Position | `GenVtx_x` | `float` | |
| Position | `GenVtx_y` | `float` | |
| Position | `GenVtx_z` | `float` | |
| Position | `GenVtx_t0` | `float` | |

## LayerCluster (14) — reco 2D per-layer clusters of RecHits; reco-level, matched to truth.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `LayerCluster_nHits` | `int32_t[]` | |
| Identity | `LayerCluster_seedDetId` | `int32_t[]` | |
| Kinematics | `LayerCluster_energy` | `float[]` | |
| Kinematics | `LayerCluster_eta` | `float[]` | |
| Kinematics | `LayerCluster_phi` | `float[]` | |
| Position | `LayerCluster_x` | `float[]` | |
| Position | `LayerCluster_y` | `float[]` | |
| Position | `LayerCluster_z` | `float[]` | |
| Matching to CaloPart | `LayerCluster_CaloPart_MatchIdx` | `int32_t[]` | |
| Matching to CaloPart | `LayerCluster_CaloPart_MatchQual` | `float[]` | |
| Matching to CaloPart | `LayerCluster_CaloPartNumMatch` | `int32_t[]` | |
| Matching to SimCluster | `LayerCluster_SimCluster_MatchIdx` | `int32_t[]` | |
| Matching to SimCluster | `LayerCluster_SimCluster_MatchQual` | `float[]` | |
| Matching to SimCluster | `LayerCluster_SimClusterNumMatch` | `int32_t[]` | |

## RecHitHGC (25) — finest-granularity reconstructed hit, with truth links to LayerCluster/SimCluster/BoundarySimCluster/MergedSimCluster. Source of the hit-overlap matching this project needs.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity/Position | `RecHitHGC_detId` | `int32_t[]` | detector-cell id |
| Identity/Position | `RecHitHGC_x` | `float[]` | |
| Identity/Position | `RecHitHGC_y` | `float[]` | |
| Identity/Position | `RecHitHGC_z` | `float[]` | |
| Identity/Position | `RecHitHGC_hitr` | `float[]` | |
| Kinematics | `RecHitHGC_energy` | `float[]` | deposited energy — basis for Jaccard/shared-energy matching |
| Kinematics | `RecHitHGC_time` | `float[]` | |
| Matching to LayerCluster | `RecHitHGC_LayerCluster_MatchIdx` | `int32_t[]` | |
| Matching to LayerCluster | `RecHitHGC_LayerCluster_MatchQual` | `float[]` | |
| Matching to LayerCluster | `RecHitHGC_LayerClusterNumMatch` | `int32_t[]` | |
| Matching to SimCluster | `RecHitHGC_SimCluster_MatchIdx` | `int32_t[]` | a hit can match multiple SimClusters |
| Matching to SimCluster | `RecHitHGC_SimCluster_MatchQual` | `float[]` | |
| Matching to SimCluster | `RecHitHGC_SimClusterNumMatch` | `int32_t[]` | |
| Matching to SimCluster | `RecHitHGC_SimClusterBestMatchIdx` | `int32_t[]` | |
| Matching to SimCluster | `RecHitHGC_SimClusterBestMatchQual` | `int32_t[]` | typed int32 despite looking like a fraction — verify values before trusting as energy fraction |
| Matching to BoundarySimCluster (*`bound_*` only*) | `RecHitHGC_BoundarySimCluster_MatchIdx` | `int32_t[]` | |
| Matching to BoundarySimCluster (*`bound_*` only*) | `RecHitHGC_BoundarySimCluster_MatchQual` | `float[]` | |
| Matching to BoundarySimCluster (*`bound_*` only*) | `RecHitHGC_BoundarySimClusterNumMatch` | `int32_t[]` | |
| Matching to BoundarySimCluster (*`bound_*` only*) | `RecHitHGC_BoundarySimClusterBestMatchIdx` | `int32_t[]` | |
| Matching to BoundarySimCluster (*`bound_*` only*) | `RecHitHGC_BoundarySimClusterBestMatchQual` | `int32_t[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimCluster_MatchIdx` | `int32_t[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimCluster_MatchQual` | `float[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimClusterNumMatch` | `int32_t[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimClusterBestMatchIdx` | `int32_t[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimClusterBestMatchQual` | `int32_t[]` | same int32 caveat as above |

## BoundarySimCluster (34) — *`bound_*` only.* The boundary collection, the central object of this project and the default (`--table boundary`). Paper concept: SimShower. One SimCluster per track that crossed the boundary with no ancestor that crossed; a descendant's hits are folded into it. Same fields as `SimCluster_*`.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Flags | `BoundarySimCluster_allHitsHGCAL` | `bool[]` |  |
| Flags | `BoundarySimCluster_crossedBoundary` | `bool[]` |  |
| Flags | `BoundarySimCluster_hasHGCALHit` | `bool[]` |  |
| Flags | `BoundarySimCluster_isTrainable` | `bool[]` |  |
| Flags | `BoundarySimCluster_onHGCFrontFace` | `bool[]` |  |
| Identity | `BoundarySimCluster_bunchCrossing` | `int32_t[]` |  |
| Identity | `BoundarySimCluster_charge` | `int32_t[]` |  |
| Identity | `BoundarySimCluster_eventId` | `int32_t[]` |  |
| Identity | `BoundarySimCluster_pdgId` | `int32_t[]` |  |
| Identity | `BoundarySimCluster_rawEventId` | `int32_t[]` |  |
| Identity | `BoundarySimCluster_trackId` | `int32_t[]` | Geant4 track id: a track stored in both files of a pair has the same id in both (`classify_decompositions.id_agreement` checks this on every run); ids restart per event |
| Identity | `BoundarySimCluster_CaloPartIdx` | `int32_t[]` | **SimCluster → CaloParticle link** |
| Boundary | `BoundarySimCluster_trackIdAtBoundary` | `int32_t[]` | the track that crossed and opened the SimCluster; compared across configs |
| Boundary | `BoundarySimCluster_boundaryEnergy` | `float[]` | energy at the crossing, rest mass included (`boundary_geometry.kinetic_energy`) |
| Boundary | `BoundarySimCluster_boundaryEnergyNoMu` | `float[]` |  |
| Boundary | `BoundarySimCluster_boundaryP4` | `float[]` |  |
| Boundary | `BoundarySimCluster_boundaryPmag` | `float[]` |  |
| Boundary | `BoundarySimCluster_boundaryPt` | `float[]` |  |
| Boundary | `BoundarySimCluster_impactPoint_eta` | `float[]` | the momentum direction, not the position: docs/data_handling.md |
| Boundary | `BoundarySimCluster_impactPoint_phi` | `float[]` |  |
| Boundary | `BoundarySimCluster_impactPoint_t` | `float[]` |  |
| Boundary | `BoundarySimCluster_impactPoint_x` | `float[]` |  |
| Boundary | `BoundarySimCluster_impactPoint_y` | `float[]` |  |
| Boundary | `BoundarySimCluster_impactPoint_z` | `float[]` | most crossings on each config's front face: Standard \|z\| = 293.5 cm, FineCalo 321.05 cm |
| Kinematics | `BoundarySimCluster_eta` | `float[]` |  |
| Kinematics | `BoundarySimCluster_mass` | `float[]` | **unusable**, as `SimCluster_mass` |
| Kinematics | `BoundarySimCluster_phi` | `float[]` |  |
| Kinematics | `BoundarySimCluster_pt` | `float[]` |  |
| Position | `BoundarySimCluster_lastPos_x` | `float[]` |  |
| Position | `BoundarySimCluster_lastPos_y` | `float[]` |  |
| Position | `BoundarySimCluster_lastPos_z` | `float[]` |  |
| Hit stats | `BoundarySimCluster_nHits` | `int32_t[]` |  |
| Hit stats | `BoundarySimCluster_sumHitEnergy` | `float[]` | hit energy of this SimCluster **and its descendants**: never sum it over SimClusters |
| Hit stats | `BoundarySimCluster_recEnergy` | `float[]` | reconstructed hit energy assigned to this SimCluster alone; summable |

## SimCluster (34) — the legacy collection: one SimCluster per track that crossed the boundary, nested ones included (`--table legacy`).

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Flags | `SimCluster_allHitsHGCAL` | `bool[]` | |
| Flags | `SimCluster_crossedBoundary` | `bool[]` | |
| Flags | `SimCluster_hasHGCALHit` | `bool[]` | |
| Flags | `SimCluster_isTrainable` | `bool[]` | |
| Flags | `SimCluster_onHGCFrontFace` | `bool[]` | ~100% true in **both** configs despite their different boundary planes — non-discriminating |
| Identity | `SimCluster_pdgId` | `int32_t[]` | |
| Identity | `SimCluster_charge` | `int32_t[]` | |
| Identity | `SimCluster_trackId` | `int32_t[]` | Geant4 track id: a track stored in both files of a pair has the same id in both (`classify_decompositions.id_agreement` checks this on every run); ids restart per event |
| Identity | `SimCluster_eventId` | `int32_t[]` | |
| Identity | `SimCluster_rawEventId` | `int32_t[]` | |
| Identity | `SimCluster_bunchCrossing` | `int32_t[]` | |
| Identity | `SimCluster_CaloPartIdx` | `int32_t[]` | **SimCluster → CaloParticle link** |
| Boundary | `SimCluster_boundaryEnergy` | `float[]` | energy at the crossing, rest mass included (`boundary_geometry.kinetic_energy`) |
| Boundary | `SimCluster_boundaryEnergyNoMu` | `float[]` | |
| Boundary | `SimCluster_boundaryP4` | `float[]` | |
| Boundary | `SimCluster_boundaryPmag` | `float[]` | |
| Boundary | `SimCluster_boundaryPt` | `float[]` | |
| Boundary | `SimCluster_trackIdAtBoundary` | `int32_t[]` | |
| Boundary | `SimCluster_impactPoint_eta` | `float[]` | the momentum direction, not the position: docs/data_handling.md |
| Boundary | `SimCluster_impactPoint_phi` | `float[]` | |
| Boundary | `SimCluster_impactPoint_t` | `float[]` | |
| Boundary | `SimCluster_impactPoint_x` | `float[]` | |
| Boundary | `SimCluster_impactPoint_y` | `float[]` | |
| Boundary | `SimCluster_impactPoint_z` | `float[]` | **two different planes per config** — fc0 \|z\|=293.5 cm, fc1 \|z\|=321.0 cm |
| Kinematics | `SimCluster_eta` | `float[]` | |
| Kinematics | `SimCluster_phi` | `float[]` | |
| Kinematics | `SimCluster_pt` | `float[]` | |
| Kinematics | `SimCluster_mass` | `float[]` | **unusable**: not a rest mass (it takes negative values); compute sqrt(E² − p²) or use `pdgId` |
| Position | `SimCluster_lastPos_x` | `float[]` | |
| Position | `SimCluster_lastPos_y` | `float[]` | |
| Position | `SimCluster_lastPos_z` | `float[]` | |
| Hit stats | `SimCluster_nHits` | `int32_t[]` | |
| Hit stats | `SimCluster_sumHitEnergy` | `float[]` | hit energy of this SimCluster **and its descendants**: never sum it over SimClusters |
| Hit stats | `SimCluster_recEnergy` | `float[]` | reconstructed hit energy assigned to this SimCluster alone; summable |

## SimTrack (19) — G4 tracks; `SimClusterIdx` / `BoundarySimClusterIdx` link a track to the SimCluster it opened.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Flags | `SimTrack_crossedBoundary` | `bool[]` | |
| Identity | `SimTrack_pdgId` | `int32_t[]` | |
| Identity | `SimTrack_charge` | `int32_t[]` | |
| Identity | `SimTrack_trackId` | `int32_t[]` | |
| Identity | `SimTrack_trackIdAtBoundary` | `int32_t[]` | |
| Identity | `SimTrack_SimClusterIdx` | `int32_t[]` | SimTrack → SimCluster (legacy) link |
| Identity | `SimTrack_BoundarySimClusterIdx` | `int32_t[]` | SimTrack → BoundarySimCluster link (*`bound_*` only*) |
| Boundary | `SimTrack_boundaryMomentum_eta` | `float[]` | momentum **at the boundary record**, not at production |
| Boundary | `SimTrack_boundaryMomentum_phi` | `float[]` | |
| Boundary | `SimTrack_boundaryMomentum_pt` | `float[]` | |
| Boundary | `SimTrack_boundaryPos_x` | `float[]` | |
| Boundary | `SimTrack_boundaryPos_y` | `float[]` | |
| Boundary | `SimTrack_boundaryPos_z` | `float[]` | |
| Kinematics | `SimTrack_eta` | `float[]` | momentum **at production** — a different vector from `boundaryMomentum_*` |
| Kinematics | `SimTrack_phi` | `float[]` | |
| Kinematics | `SimTrack_pt` | `float[]` | |
| Position | `SimTrack_lastPos_x` | `float[]` | |
| Position | `SimTrack_lastPos_y` | `float[]` | |
| Position | `SimTrack_lastPos_z` | `float[]` | |

## Event-level (31) — per-event scalars, not per-object arrays.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Event identifiers | `run` | `uint32_t` | does **not** identify matching events across fc0/fc1 — see docs/data_handling.md |
| Event identifiers | `luminosityBlock` | `uint32_t` | |
| Event identifiers | `event` | `uint64_t` | |
| Event identifiers | `bunchCrossing` | `uint32_t` | |
| Event identifiers | `orbitNumber` | `uint32_t` | |
| Collection sizes | `nCaloPart` | `int32_t` | |
| Collection sizes | `nGenPart` | `int32_t` | |
| Collection sizes | `nTrack` | `int32_t` | |
| Collection sizes | `nSimHitHGCEE` | `int32_t` | |
| Collection sizes | `nSimHitHGCHEB` | `int32_t` | |
| Collection sizes | `nSimHitHGCHEF` | `int32_t` | |
| Collection sizes | `nRecHitHGC` | `int32_t` | |
| Collection sizes | `nLayerCluster` | `int32_t` | |
| Collection sizes | `nMergedSimCluster` | `int32_t` | |
| Collection sizes | `nPFCand` | `int32_t` | |
| Collection sizes | `nPFTICLCand` | `int32_t` | |
| Collection sizes | `nSimCluster` | `int32_t` | |
| Collection sizes | `nBoundarySimCluster` | `int32_t` | *`bound_*` only* |
| Collection sizes | `nSimTrack` | `int32_t` | |
| Collection sizes | `nTrackConv` | `int32_t` | |
| Collection sizes | `nTrackDisp` | `int32_t` | |
| Collection sizes | `nTrackingPart` | `int32_t` | |
| Cross-collection relationship counts | `nLayerCluster_CaloPart` | `int32_t` | |
| Cross-collection relationship counts | `nLayerCluster_SimCluster` | `int32_t` | |
| Cross-collection relationship counts | `nMergedSimCluster_SimCluster` | `int32_t` | |
| Cross-collection relationship counts | `nRecHitHGC_LayerCluster` | `int32_t` | |
| Cross-collection relationship counts | `nRecHitHGC_MergedSimCluster` | `int32_t` | |
| Cross-collection relationship counts | `nRecHitHGC_SimCluster` | `int32_t` | |
| Cross-collection relationship counts | `nRecHitHGC_BoundarySimCluster` | `int32_t` | *`bound_*` only* |
| Cross-collection relationship counts | `nTrack_TrackingPart` | `int32_t` | |
| Cross-collection relationship counts | `nTrackingPart_Track` | `int32_t` | |

---

**Total: 414 branches** (17+9+4+14+36+9+9+25+34+34+7+7+7+19+84+19+16+8+25+31); 288 in the files without the boundary collection.
## Reference only — collections not needed for the boundary-SimCluster question

Present in the files and documented for completeness, but none of these feeds
the CaloParticle → boundary-SimCluster comparison:

- `MergedSimCluster_*` — post-merge, downstream of the decomposition being
  validated.
- `PFCand_*`, `PFTICLCand_*`, `Track*`, `TrackingPart_*` — reco level.
  `TrackConv_*` is the one plausible future exception: it refits tracks under
  the photon-conversion hypothesis, which is the topology these comparisons are
  about.
- `TICLCand_*`, `ticlTracksters*`, `Candidate2Tracksters*` — TICL reconstruction output
  (`bound_*` only).
- `SimHitHGC{EE,HEB,HEF}_*` — raw G4 deposits one level below `RecHitHGC`,
  needed only if `RecHitHGC`'s SimCluster match proves too coarse.

## MergedSimCluster (36) — *not needed for this project's core question.* Post-merge target object; downstream of the boundary-SimCluster decomposition being validated.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Flags | `MergedSimCluster_allHitsHGCAL` | `bool[]` | |
| Flags | `MergedSimCluster_crossedBoundary` | `bool[]` | |
| Flags | `MergedSimCluster_hasHGCALHit` | `bool[]` | |
| Flags | `MergedSimCluster_isTrainable` | `bool[]` | |
| Flags | `MergedSimCluster_onHGCFrontFace` | `bool[]` | |
| Identity | `MergedSimCluster_pdgId` | `int32_t[]` | |
| Identity | `MergedSimCluster_charge` | `int32_t[]` | |
| Identity | `MergedSimCluster_trackId` | `int32_t[]` | |
| Identity | `MergedSimCluster_eventId` | `int32_t[]` | |
| Identity | `MergedSimCluster_rawEventId` | `int32_t[]` | |
| Identity | `MergedSimCluster_bunchCrossing` | `int32_t[]` | |
| Boundary | `MergedSimCluster_boundaryEnergy` | `float[]` | |
| Boundary | `MergedSimCluster_boundaryEnergyNoMu` | `float[]` | |
| Boundary | `MergedSimCluster_boundaryP4` | `float[]` | |
| Boundary | `MergedSimCluster_boundaryPmag` | `float[]` | |
| Boundary | `MergedSimCluster_boundaryPt` | `float[]` | |
| Boundary | `MergedSimCluster_trackIdAtBoundary` | `int32_t[]` | |
| Boundary | `MergedSimCluster_impactPoint_eta` | `float[]` | |
| Boundary | `MergedSimCluster_impactPoint_phi` | `float[]` | |
| Boundary | `MergedSimCluster_impactPoint_t` | `float[]` | |
| Boundary | `MergedSimCluster_impactPoint_x` | `float[]` | |
| Boundary | `MergedSimCluster_impactPoint_y` | `float[]` | |
| Boundary | `MergedSimCluster_impactPoint_z` | `float[]` | |
| Kinematics | `MergedSimCluster_eta` | `float[]` | |
| Kinematics | `MergedSimCluster_phi` | `float[]` | |
| Kinematics | `MergedSimCluster_pt` | `float[]` | |
| Kinematics | `MergedSimCluster_mass` | `float[]` | |
| Position | `MergedSimCluster_lastPos_x` | `float[]` | |
| Position | `MergedSimCluster_lastPos_y` | `float[]` | |
| Position | `MergedSimCluster_lastPos_z` | `float[]` | |
| Hit stats | `MergedSimCluster_nHits` | `int32_t[]` | |
| Hit stats | `MergedSimCluster_sumHitEnergy` | `float[]` | |
| Hit stats | `MergedSimCluster_recEnergy` | `float[]` | |
| Matching to SimCluster | `MergedSimCluster_SimCluster_MatchIdx` | `int32_t[]` | |
| Matching to SimCluster | `MergedSimCluster_SimCluster_MatchQual` | `float[]` | |
| Matching to SimCluster | `MergedSimCluster_SimClusterNumMatch` | `int32_t[]` | |

## PFCand (9) — *reco-level, not needed for this project's core question.* Final particle-flow candidates.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `PFCand_pdgId` | `int32_t[]` | |
| Identity | `PFCand_charge` | `int32_t[]` | |
| Position | `PFCand_Vtx_x` | `float[]` | |
| Position | `PFCand_Vtx_y` | `float[]` | |
| Position | `PFCand_Vtx_z` | `float[]` | |
| Kinematics | `PFCand_eta` | `float[]` | |
| Kinematics | `PFCand_phi` | `float[]` | |
| Kinematics | `PFCand_pt` | `float[]` | |
| Kinematics | `PFCand_mass` | `float[]` | |

## PFTICLCand (9) — *reco-level, not needed for this project's core question.* TICL-based particle-flow candidates. Same field structure as PFCand.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `PFTICLCand_pdgId` | `int32_t[]` | |
| Identity | `PFTICLCand_charge` | `int32_t[]` | |
| Position | `PFTICLCand_Vtx_x` | `float[]` | |
| Position | `PFTICLCand_Vtx_y` | `float[]` | |
| Position | `PFTICLCand_Vtx_z` | `float[]` | |
| Kinematics | `PFTICLCand_eta` | `float[]` | |
| Kinematics | `PFTICLCand_phi` | `float[]` | |
| Kinematics | `PFTICLCand_pt` | `float[]` | |
| Kinematics | `PFTICLCand_mass` | `float[]` | |

## SimHitHGCEE / SimHitHGCHEB / SimHitHGCHEF (7 each) — *not needed unless RecHitHGC's SimCluster match proves too coarse.* Raw G4 energy deposits, pre-reconstruction, per calorimeter subdetector (EE, HE-back, HE-front).

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `SimHitHGC{EE,HEB,HEF}_detId` | `int32_t[]` | |
| Identity | `SimHitHGC{EE,HEB,HEF}_trackId` | `int32_t[]` | |
| Kinematics | `SimHitHGC{EE,HEB,HEF}_energy` | `float[]` | |
| Position | `SimHitHGC{EE,HEB,HEF}_x` | `float[]` | |
| Position | `SimHitHGC{EE,HEB,HEF}_y` | `float[]` | |
| Position | `SimHitHGC{EE,HEB,HEF}_z` | `float[]` | |
| Position | `SimHitHGC{EE,HEB,HEF}_hitr` | `float[]` | |

## Track (19) — *reco-level, not needed for this project's core question.* Reconstructed charged tracks, matched to TrackingPart truth.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `Track_charge` | `int32_t[]` | |
| Identity | `Track_numberOfLostHits` | `int32_t[]` | |
| Identity | `Track_numberOfValidHits` | `int32_t[]` | |
| Position | `Track_DecayVtx_t` | `float[]` | |
| Position | `Track_DecayVtx_x` | `float[]` | |
| Position | `Track_DecayVtx_y` | `float[]` | |
| Position | `Track_DecayVtx_z` | `float[]` | |
| Position | `Track_Vtx_t` | `float[]` | |
| Position | `Track_Vtx_x` | `float[]` | |
| Position | `Track_Vtx_y` | `float[]` | |
| Position | `Track_Vtx_z` | `float[]` | |
| Kinematics | `Track_eta` | `float[]` | |
| Kinematics | `Track_normChiSq` | `float[]` | |
| Kinematics | `Track_p` | `float[]` | |
| Kinematics | `Track_phi` | `float[]` | |
| Kinematics | `Track_pt` | `float[]` | |
| Matching to TrackingPart | `Track_TrackingPart_MatchIdx` | `int32_t[]` | |
| Matching to TrackingPart | `Track_TrackingPart_MatchQual` | `float[]` | |
| Matching to TrackingPart | `Track_TrackingPartNumMatch` | `int32_t[]` | |

## TrackConv (16) — *reco-level, not needed for this project's core question.* Tracks refit under the photon-conversion hypothesis — could help identify photon-conversion topologies.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `TrackConv_charge` | `int32_t[]` | |
| Identity | `TrackConv_numberOfLostHits` | `int32_t[]` | |
| Identity | `TrackConv_numberOfValidHits` | `int32_t[]` | |
| Position | `TrackConv_DecayVtx_t` | `float[]` | |
| Position | `TrackConv_DecayVtx_x` | `float[]` | |
| Position | `TrackConv_DecayVtx_y` | `float[]` | |
| Position | `TrackConv_DecayVtx_z` | `float[]` | |
| Position | `TrackConv_Vtx_t` | `float[]` | |
| Position | `TrackConv_Vtx_x` | `float[]` | |
| Position | `TrackConv_Vtx_y` | `float[]` | |
| Position | `TrackConv_Vtx_z` | `float[]` | |
| Kinematics | `TrackConv_eta` | `float[]` | |
| Kinematics | `TrackConv_normChiSq` | `float[]` | |
| Kinematics | `TrackConv_p` | `float[]` | |
| Kinematics | `TrackConv_phi` | `float[]` | |
| Kinematics | `TrackConv_pt` | `float[]` | |

## TrackDisp (8) — *reco-level, not needed for this project's core question.* Tracks refit under a displaced-vertex hypothesis.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `TrackDisp_charge` | `int32_t[]` | |
| Position | `TrackDisp_Vtx_t` | `float[]` | |
| Position | `TrackDisp_Vtx_x` | `float[]` | |
| Position | `TrackDisp_Vtx_y` | `float[]` | |
| Position | `TrackDisp_Vtx_z` | `float[]` | |
| Kinematics | `TrackDisp_eta` | `float[]` | |
| Kinematics | `TrackDisp_phi` | `float[]` | |
| Kinematics | `TrackDisp_pt` | `float[]` | |

## TrackingPart (25) — *reco-level, not needed for this project's core question.* GEN/SIM-level truth object matched to reconstructed Tracks.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Identity | `TrackingPart_pdgId` | `int32_t[]` | |
| Identity | `TrackingPart_charge` | `int32_t[]` | |
| Identity | `TrackingPart_trackId` | `int32_t[]` | |
| Identity | `TrackingPart_GenPartIdx` | `int32_t[]` | |
| Identity | `TrackingPart_eventId` | `int32_t[]` | |
| Identity | `TrackingPart_rawEventId` | `int32_t[]` | |
| Identity | `TrackingPart_bunchCrossing` | `int32_t[]` | |
| Identity | `TrackingPart_nDecayVtx` | `int32_t[]` | |
| Identity | `TrackingPart_nGenPart` | `int32_t[]` | |
| Identity | `TrackingPart_nSimTrack` | `int32_t[]` | |
| Position | `TrackingPart_DecayVtx_t` | `float[]` | |
| Position | `TrackingPart_DecayVtx_x` | `float[]` | |
| Position | `TrackingPart_DecayVtx_y` | `float[]` | |
| Position | `TrackingPart_DecayVtx_z` | `float[]` | |
| Position | `TrackingPart_Vtx_t` | `float[]` | |
| Position | `TrackingPart_Vtx_x` | `float[]` | |
| Position | `TrackingPart_Vtx_y` | `float[]` | |
| Position | `TrackingPart_Vtx_z` | `float[]` | |
| Kinematics | `TrackingPart_eta` | `float[]` | |
| Kinematics | `TrackingPart_mass` | `float[]` | |
| Kinematics | `TrackingPart_phi` | `float[]` | |
| Kinematics | `TrackingPart_pt` | `float[]` | |
| Matching to Track | `TrackingPart_Track_MatchIdx` | `int32_t[]` | |
| Matching to Track | `TrackingPart_Track_MatchQual` | `float[]` | |
| Matching to Track | `TrackingPart_TrackNumMatch` | `int32_t[]` | |

## TICL (84) — *`bound_*` only; reco-level, not needed for this project's core question.* TICL candidates and the tracksters of four TICL steps (CLUE3DHigh, Recovery, Links, SuperclusteringDNN), with their counters and vertex lists.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Collection sizes | `nCandidate2Tracksters` | `int32_t` | |
| Candidate2Tracksters | `Candidate2Tracksters_nCandidate2TrackstersIndices` | `uint16_t[]` | |
| Collection sizes | `nCandidate2TrackstersIndices` | `int32_t` | |
| Candidate2TrackstersIndices | `Candidate2TrackstersIndices_tracksterIndex` | `uint32_t[]` | |
| Collection sizes | `nTICLCand` | `int32_t` | |
| TICLCand | `TICLCand_charge` | `int32_t[]` | |
| TICLCand | `TICLCand_pdgID` | `int32_t[]` | |
| TICLCand | `TICLCand_energy` | `float[]` | |
| TICLCand | `TICLCand_eta` | `float[]` | |
| TICLCand | `TICLCand_mass` | `float[]` | |
| TICLCand | `TICLCand_p` | `float[]` | |
| TICLCand | `TICLCand_phi` | `float[]` | |
| TICLCand | `TICLCand_pt` | `float[]` | |
| TICLCand | `TICLCand_raw_energy` | `float[]` | |
| TICLCand | `TICLCand_time` | `float[]` | |
| TICLCand | `TICLCand_timeError` | `float[]` | |
| Collection sizes | `nticlTracksterLinksSuperclusteringDNN` | `int32_t` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_nticlTracksterLinksSuperclusteringDNNvertices` | `uint16_t[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_oticlTracksterLinksSuperclusteringDNNvertices` | `uint16_t[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_barycenter_eta` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_barycenter_phi` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_barycenter_x` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_barycenter_y` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_barycenter_z` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_raw_em_energy` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_raw_energy` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_raw_pt` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_regressed_energy` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_time` | `float[]` | |
| ticlTracksterLinksSuperclusteringDNN | `ticlTracksterLinksSuperclusteringDNN_timeError` | `float[]` | |
| Collection sizes | `nticlTracksterLinksSuperclusteringDNNvertices` | `int32_t` | |
| ticlTracksterLinksSuperclusteringDNNvertices | `ticlTracksterLinksSuperclusteringDNNvertices_vertices` | `uint32_t[]` | |
| ticlTracksterLinksSuperclusteringDNNvertices | `ticlTracksterLinksSuperclusteringDNNvertices_vertex_mult` | `float[]` | |
| Collection sizes | `nticlTracksterLinks` | `int32_t` | |
| ticlTracksterLinks | `ticlTracksterLinks_nticlTracksterLinksvertices` | `uint16_t[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_oticlTracksterLinksvertices` | `uint16_t[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_barycenter_eta` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_barycenter_phi` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_barycenter_x` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_barycenter_y` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_barycenter_z` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_raw_em_energy` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_raw_energy` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_raw_pt` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_regressed_energy` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_time` | `float[]` | |
| ticlTracksterLinks | `ticlTracksterLinks_timeError` | `float[]` | |
| Collection sizes | `nticlTracksterLinksvertices` | `int32_t` | |
| ticlTracksterLinksvertices | `ticlTracksterLinksvertices_vertices` | `uint32_t[]` | |
| ticlTracksterLinksvertices | `ticlTracksterLinksvertices_vertex_mult` | `float[]` | |
| Collection sizes | `nticlTrackstersCLUE3DHigh` | `int32_t` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_nticlTrackstersCLUE3DHighvertices` | `uint16_t[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_oticlTrackstersCLUE3DHighvertices` | `uint16_t[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_barycenter_eta` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_barycenter_phi` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_barycenter_x` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_barycenter_y` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_barycenter_z` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_raw_em_energy` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_raw_energy` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_raw_pt` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_regressed_energy` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_time` | `float[]` | |
| ticlTrackstersCLUE3DHigh | `ticlTrackstersCLUE3DHigh_timeError` | `float[]` | |
| Collection sizes | `nticlTrackstersCLUE3DHighvertices` | `int32_t` | |
| ticlTrackstersCLUE3DHighvertices | `ticlTrackstersCLUE3DHighvertices_vertices` | `uint32_t[]` | |
| ticlTrackstersCLUE3DHighvertices | `ticlTrackstersCLUE3DHighvertices_vertex_mult` | `float[]` | |
| Collection sizes | `nticlTrackstersRecovery` | `int32_t` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_nticlTrackstersRecoveryvertices` | `uint16_t[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_oticlTrackstersRecoveryvertices` | `uint16_t[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_barycenter_eta` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_barycenter_phi` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_barycenter_x` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_barycenter_y` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_barycenter_z` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_raw_em_energy` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_raw_energy` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_raw_pt` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_regressed_energy` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_time` | `float[]` | |
| ticlTrackstersRecovery | `ticlTrackstersRecovery_timeError` | `float[]` | |
| Collection sizes | `nticlTrackstersRecoveryvertices` | `int32_t` | |
| ticlTrackstersRecoveryvertices | `ticlTrackstersRecoveryvertices_vertices` | `uint32_t[]` | |
| ticlTrackstersRecoveryvertices | `ticlTrackstersRecoveryvertices_vertex_mult` | `float[]` | |
