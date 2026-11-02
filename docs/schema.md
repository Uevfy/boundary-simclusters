# ROOT/nanoML schema reference

All 288 branches of the `Events` tree, grouped by collection and by relation
within each collection. Confirmed identical (name, type, count) across all
original local test files (`testNanoML.root`, `FinetestNanoML.root` and a
second 200 GeV e⁻ pair since deleted — 1000 events each).
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
| Identity | `CaloPart_trackId` | `int32_t[]` | G4 track id — not stable across fc0/fc1 |
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

## RecHitHGC (20) — finest-granularity reconstructed hit, with truth links to LayerCluster/SimCluster/MergedSimCluster. Source of the hit-overlap matching this project needs.

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
| Matching to MergedSimCluster | `RecHitHGC_MergedSimCluster_MatchIdx` | `int32_t[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimCluster_MatchQual` | `float[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimClusterNumMatch` | `int32_t[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimClusterBestMatchIdx` | `int32_t[]` | |
| Matching to MergedSimCluster | `RecHitHGC_MergedSimClusterBestMatchQual` | `int32_t[]` | same int32 caveat as above |

## SimCluster (34) — boundary SimCluster, the central object of this project. Paper concept: SimShower.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Flags | `SimCluster_allHitsHGCAL` | `bool[]` | |
| Flags | `SimCluster_crossedBoundary` | `bool[]` | |
| Flags | `SimCluster_hasHGCALHit` | `bool[]` | |
| Flags | `SimCluster_isTrainable` | `bool[]` | |
| Flags | `SimCluster_onHGCFrontFace` | `bool[]` | ~100% true in **both** configs despite their different boundary planes — non-discriminating |
| Identity | `SimCluster_pdgId` | `int32_t[]` | |
| Identity | `SimCluster_charge` | `int32_t[]` | |
| Identity | `SimCluster_trackId` | `int32_t[]` | G4 track id — not stable across fc0/fc1, do not match on this |
| Identity | `SimCluster_eventId` | `int32_t[]` | |
| Identity | `SimCluster_rawEventId` | `int32_t[]` | |
| Identity | `SimCluster_bunchCrossing` | `int32_t[]` | |
| Identity | `SimCluster_CaloPartIdx` | `int32_t[]` | **SimCluster → CaloParticle link** |
| Boundary | `SimCluster_boundaryEnergy` | `float[]` | |
| Boundary | `SimCluster_boundaryEnergyNoMu` | `float[]` | |
| Boundary | `SimCluster_boundaryP4` | `float[]` | |
| Boundary | `SimCluster_boundaryPmag` | `float[]` | |
| Boundary | `SimCluster_boundaryPt` | `float[]` | |
| Boundary | `SimCluster_trackIdAtBoundary` | `int32_t[]` | |
| Boundary | `SimCluster_impactPoint_eta` | `float[]` | |
| Boundary | `SimCluster_impactPoint_phi` | `float[]` | |
| Boundary | `SimCluster_impactPoint_t` | `float[]` | |
| Boundary | `SimCluster_impactPoint_x` | `float[]` | |
| Boundary | `SimCluster_impactPoint_y` | `float[]` | |
| Boundary | `SimCluster_impactPoint_z` | `float[]` | **two different planes per config** — fc0 \|z\|=293.5 cm, fc1 \|z\|=321.0 cm |
| Kinematics | `SimCluster_eta` | `float[]` | |
| Kinematics | `SimCluster_phi` | `float[]` | |
| Kinematics | `SimCluster_pt` | `float[]` | |
| Kinematics | `SimCluster_mass` | `float[]` | **unusable** — see docs/data_handling.md |
| Position | `SimCluster_lastPos_x` | `float[]` | |
| Position | `SimCluster_lastPos_y` | `float[]` | |
| Position | `SimCluster_lastPos_z` | `float[]` | |
| Hit stats | `SimCluster_nHits` | `int32_t[]` | |
| Hit stats | `SimCluster_sumHitEnergy` | `float[]` | one of three energy definitions — see docs/data_handling.md |
| Hit stats | `SimCluster_recEnergy` | `float[]` | |

## SimTrack (18) — G4 track segments; `SimClusterIdx` links a track to the boundary SimCluster it belongs to.

| Subgroup | Field | Type | Note |
|---|---|---|---|
| Flags | `SimTrack_crossedBoundary` | `bool[]` | |
| Identity | `SimTrack_pdgId` | `int32_t[]` | |
| Identity | `SimTrack_charge` | `int32_t[]` | |
| Identity | `SimTrack_trackId` | `int32_t[]` | |
| Identity | `SimTrack_trackIdAtBoundary` | `int32_t[]` | |
| Identity | `SimTrack_SimClusterIdx` | `int32_t[]` | SimTrack → SimCluster link |
| Boundary | `SimTrack_boundaryMomentum_eta` | `float[]` | momentum **at the boundary record**, not at production |
| Boundary | `SimTrack_boundaryMomentum_phi` | `float[]` | see docs/data_handling.md |
| Boundary | `SimTrack_boundaryMomentum_pt` | `float[]` | see docs/data_handling.md |
| Boundary | `SimTrack_boundaryPos_x` | `float[]` | |
| Boundary | `SimTrack_boundaryPos_y` | `float[]` | |
| Boundary | `SimTrack_boundaryPos_z` | `float[]` | |
| Kinematics | `SimTrack_eta` | `float[]` | momentum **at production** — a different vector from `boundaryMomentum_*` |
| Kinematics | `SimTrack_phi` | `float[]` | see docs/data_handling.md |
| Kinematics | `SimTrack_pt` | `float[]` | see docs/data_handling.md |
| Position | `SimTrack_lastPos_x` | `float[]` | |
| Position | `SimTrack_lastPos_y` | `float[]` | |
| Position | `SimTrack_lastPos_z` | `float[]` | |

## Event-level (29) — per-event scalars, not per-object arrays.

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
| Cross-collection relationship counts | `nTrack_TrackingPart` | `int32_t` | |
| Cross-collection relationship counts | `nTrackingPart_Track` | `int32_t` | |

---

**Total: 288 branches** (17+9+4+14+36+9+9+20+34+7+7+7+18+19+16+8+25+29).
## Reference only — collections not needed for the boundary-SimCluster question

Present in the files and documented for completeness, but none of these feeds
the CaloParticle → boundary-SimCluster comparison:

- `MergedSimCluster_*` — post-merge, downstream of the decomposition being
  validated.
- `PFCand_*`, `PFTICLCand_*`, `Track*`, `TrackingPart_*` — reco level.
  `TrackConv_*` is the one plausible future exception: it refits tracks under
  the photon-conversion hypothesis, which is the topology these comparisons are
  about.
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

