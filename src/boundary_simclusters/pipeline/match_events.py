"""Diagnose and rectify fc0/fc1 event misalignment.

fc0 and fc1 hold the same physical events, but not in the same row order —
and the standard `run`/`luminosityBlock`/`event` branches cannot be used to
join them (see docs/data_handling.md). This step:

1. Diagnoses the misalignment: shows that naive row order and the `event`
   branch both fail to pair the correct events, validated against an
   independent field (LayerCluster count per event, which must agree for
   the same physical event but is not used to build any of the pairings).
2. Rectifies it: pairs events on their generator record -- the sorted
   (pdgId, pt, eta, phi) of every GenPart. Both configs start from the same
   GEN step, so the key is bit-identical for a true pair and needs no
   tolerance. The guns draw eta and phi per event (pt follows from the fixed
   energy), which keeps the key unique (1000/1000 on every local pair,
   docs/data_handling.md).
   Events whose key is not unique in both files are left unmatched rather
   than guessed at.
3. Cross-validates the GEN mapping against a physics fingerprint (hit count +
   total deposited energy in RecHitHGC), which is what this script paired on
   at first. If a future gun fixes every generator variable, the GEN key
   collapses and the fingerprint is the fallback.
4. Pairs the CaloParticles within every matched event on exact (pdgId, eta,
   phi). Their stored order is not the same in both files
   (docs/data_handling.md), so no script may assume index `cp` is the same
   particle on both sides. Raises unless every CaloParticle pairs 1:1.
5. Reports events left unmatched rather than forcing or silently dropping
   them.
6. Writes the validated mappings to data/processed/event_matching_<tag>.parquet
   and calopart_matching_<tag>.parquet, each carrying the names of the two
   files it was built from (`read_mapping` refuses any other pair).

The match step of scripts/analyze.py (`--steps match`), `run(sample)` for the
tests.
"""

from boundary_simclusters.matching import build_event_mapping, diagnose_misalignment, save_mapping


def run(sample):
    """Pair `sample`'s events and CaloParticles, write the mappings, and return
    (event mapping, CaloParticle mapping) as written."""
    fc0_path, fc1_path = sample.files["fc0"], sample.files["fc1"]
    t0, t1 = sample.tree("fc0"), sample.tree("fc1")
    print(f"fc0: {fc0_path.name} ({t0.num_entries} events)")
    print(f"fc1: {fc1_path.name} ({t1.num_entries} events)")
    print(f"tag: {sample.tag}\n")

    diagnose_misalignment(t0, t1)
    mapping, caloparts, unmatched_fc0, unmatched_fc1 = build_event_mapping(t0, t1)
    return save_mapping(mapping, caloparts, unmatched_fc0, unmatched_fc1, sample.tag,
                        fc0_path, fc1_path)
