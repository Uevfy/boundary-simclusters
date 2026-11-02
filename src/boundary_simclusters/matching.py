"""Pairing an fc0/fc1 file pair: events on their generator record, then the
CaloParticles within each event. The library behind `pipeline.match_events`.

fc0 and fc1 hold the same physical events, but not in the same row order, and
the `run`/`luminosityBlock`/`event` branches cannot join them (docs/data_handling.md). Events are paired on their generator record -- the sorted (pdgId, pt, eta,
phi) of every GenPart, bit-identical for a true pair -- cross-checked
against a RecHit fingerprint, and the CaloParticles of each matched event on
exact (pdgId, eta, phi), since their stored order differs between the files in
some events. Repeated keys are left unmatched rather than guessed at.

`save_mapping` writes event_matching_<tag> / calopart_matching_<tag> /
unmatched_events_<tag> to data/processed/, each carrying the two file names it
was built from; `read_mapping` refuses a mapping for any other pair.
"""

from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from boundary_simclusters import output_paths



def _lengths(tree, branch):
    return np.array([len(x) for x in tree[branch].array(library="np")])


def _totals(tree, branch):
    return np.array([float(x.sum()) for x in tree[branch].array(library="np")])


def greedy_match(keys0, keys1):
    """1:1 greedy match on an exact key. Returns fc1 row per fc0 row, -1 if none."""
    pool = defaultdict(list)
    for row, key in enumerate(keys1):
        pool[key].append(row)

    out = np.full(len(keys0), -1, dtype=int)
    for row, key in enumerate(keys0):
        if pool.get(key):
            out[row] = pool[key].pop(0)
    return out


# Same-event totals can differ by ~0.007 GeV out of ~820 between the two configs,
# purely from float32 summation order. Exact equality on the energy would discard
# valid pairs, so it is matched within a tolerance well above that noise and well
# below the spacing between genuinely different events.
ENERGY_TOL = 0.02


def hit_content(tree):
    return _lengths(tree, "RecHitHGC_energy"), _totals(tree, "RecHitHGC_energy")


def match_by_hit_content(tree0, tree1, energy_tol=ENERGY_TOL):
    """Pair events on hit count (exact -- identical for true pairs) plus total
    deposited energy (within tolerance). Best candidates are consumed first, so
    a near-tie cannot steal a better pair's partner."""
    n0, e0 = hit_content(tree0)
    n1, e1 = hit_content(tree1)

    by_count = defaultdict(list)
    for row, count in enumerate(n1):
        by_count[count].append(row)

    candidates = []
    for row0, count in enumerate(n0):
        for row1 in by_count.get(count, []):
            delta = abs(e0[row0] - e1[row1])
            if delta <= energy_tol:
                candidates.append((delta, row0, row1))
    candidates.sort()

    out = np.full(len(n0), -1, dtype=int)
    used0, used1 = set(), set()
    for _, row0, row1 in candidates:
        if row0 in used0 or row1 in used1:
            continue
        out[row0] = row1
        used0.add(row0)
        used1.add(row1)
    return out


def genpart_signature(tree):
    cols = [tree[f"GenPart_{b}"].array(library="np") for b in ("pdgId", "pt", "eta", "phi")]
    return [tuple(sorted(zip(*(c.tolist() for c in event)))) for event in zip(*cols)]


def match_by_genpart(tree0, tree1):
    """Pair events on an exact GEN key, using only keys that occur once in each
    file -- a repeated key cannot say which of its events is which."""
    keys0, keys1 = genpart_signature(tree0), genpart_signature(tree1)
    count0, count1 = Counter(keys0), Counter(keys1)
    row1_of = {key: row for row, key in enumerate(keys1) if count1[key] == 1}

    out = np.full(len(keys0), -1, dtype=int)
    for row0, key in enumerate(keys0):
        if count0[key] == 1 and key in row1_of:
            out[row0] = row1_of[key]
    ambiguous = sum(n for n in count0.values() if n > 1)
    return out, ambiguous


def calopart_keys(tree):
    """Per event, the (pdgId, eta, phi) of each CaloParticle, in stored order."""
    cols = [tree[f"CaloPart_{b}"].array(library="np") for b in ("pdgId", "eta", "phi")]
    return [list(zip(*(c.tolist() for c in event))) for event in zip(*cols)]


def pair_caloparticles(keys0, keys1):
    """fc1 index of each fc0 CaloParticle in one matched event, or None.

    Paired on exact (pdgId, eta, phi): both configs start from the same GEN
    step, so a CaloParticle's kinematics are bit-identical across the pair,
    as the event key is. The stored order is not (docs/data_handling.md), and
    `CaloPart_GenPartIdx` is -1 throughout, so the kinematics are the key.
    None unless every CaloParticle finds exactly one partner.
    """
    if len(set(keys0)) != len(keys0) or sorted(keys0) != sorted(keys1):
        return None
    position = {key: cp for cp, key in enumerate(keys1)}
    return [position[key] for key in keys0]


def match_caloparticles(t0, t1, mapping):
    """One row per CaloParticle of every matched event, fc0 index -> fc1 index.

    Raises if any matched event's CaloParticles do not pair 1:1: then the event
    mapping itself is wrong, and nothing downstream can be trusted."""
    keys0, keys1 = calopart_keys(t0), calopart_keys(t1)
    rows, bad, swapped = [], [], 0
    for row0 in np.where(mapping >= 0)[0]:
        row1 = int(mapping[row0])
        partner = pair_caloparticles(keys0[row0], keys1[row1])
        if partner is None:
            bad.append((int(row0), row1))
            continue
        swapped += partner != list(range(len(partner)))
        rows += [(int(row0), cp0, row1, cp1) for cp0, cp1 in enumerate(partner)]
    if bad:
        raise RuntimeError(
            f"CaloParticles do not pair 1:1 for {len(bad)} matched event(s): {bad[:5]}... "
            "the mapping is not trustworthy for this file pair."
        )
    table = pd.DataFrame(rows, columns=["fc0_row", "fc0_cp", "fc1_row", "fc1_cp"])
    return table, swapped


# --- 1. diagnose ------------------------------------------------------------


def diagnose_misalignment(t0, t1):
    """Show that naive row order and the `event` branch fail to pair events,
    validated against LayerCluster count (independent of all three pairings)."""
    n = t0.num_entries
    lc0, lc1 = _lengths(t0, "LayerCluster_energy"), _lengths(t1, "LayerCluster_energy")
    event0 = [int(x) for x in t0["event"].array(library="np")]
    event1 = [int(x) for x in t1["event"].array(library="np")]

    pairings = {
        "naive row order": np.arange(n),
        "event branch": greedy_match(event0, event1),
        "GenPart kinematics": match_by_genpart(t0, t1)[0],
    }

    rows = {}
    for name, mapping in pairings.items():
        paired = mapping >= 0
        agree = lc0[paired] == lc1[mapping[paired]]
        rows[name] = {
            "pairs formed": int(paired.sum()),
            "LayerCluster count agrees": int(agree.sum()),
            "agreement": f"{100 * agree.mean():.1f}%",
        }

    table = pd.DataFrame(rows).T
    print("Misalignment diagnosis (validated against LayerCluster count, an")
    print("independent field not used to build any of these three pairings):\n")
    print(table.to_string())
    print()
    return table


# --- 2 & 3. rectify + cross-validate ----------------------------------------


def build_event_mapping(t0, t1):
    """Build the fc0->fc1 mapping from the generator record, cross-validate it
    against the RecHit fingerprint (a signal the GEN key never sees), and pair
    the CaloParticles within every matched event."""
    n = t0.num_entries
    mapping, ambiguous = match_by_genpart(t0, t1)

    # cross-check 1: independent matcher built from detector-level info
    hit_mapping = match_by_hit_content(t0, t1)
    both = (mapping >= 0) & (hit_mapping >= 0)
    hit_agree = int((mapping[both] == hit_mapping[both]).sum())
    if hit_agree != int(both.sum()):
        raise RuntimeError(
            f"GenPart and hit-content matchers disagree on {int(both.sum()) - hit_agree} "
            "event(s) -- the mapping is not trustworthy for this file pair."
        )

    # cross-check 2, and the CaloParticle pairing: every CaloParticle of a
    # matched event must find exactly one partner on exact kinematics
    caloparts, swapped = match_caloparticles(t0, t1, mapping)

    unmatched_fc0 = np.where(mapping < 0)[0]
    unmatched_fc1 = np.setdiff1d(np.arange(n), mapping[mapping >= 0])

    print(
        f"Rectified mapping: {int((mapping >= 0).sum())}/{n} events matched.\n"
        f"  cross-check vs RecHit fingerprint matcher   : agree on all {hit_agree} common pairs\n"
        f"  CaloParticles paired on exact kinematics    : {len(caloparts)} in all "
        f"{int((mapping >= 0).sum())} pairs\n"
        f"  events with CaloParticles in another order  : {swapped}\n"
        f"  fc0 events with a repeated GEN key          : {ambiguous}\n"
        f"  unmatched (no counterpart found)            : {len(unmatched_fc0)} fc0 / "
        f"{len(unmatched_fc1)} fc1 rows"
    )
    if len(unmatched_fc0):
        print(f"    fc0 rows: {unmatched_fc0.tolist()}")
        print(f"    fc1 rows: {unmatched_fc1.tolist()}")
    print()

    return mapping, caloparts, unmatched_fc0, unmatched_fc1


def save_mapping(mapping, caloparts, unmatched_fc0, unmatched_fc1, tag, fc0_path, fc1_path):
    out_dir = output_paths.PROCESSED
    out_dir.mkdir(parents=True, exist_ok=True)

    matched = mapping >= 0
    table = pd.DataFrame(
        {
            "fc0_row": np.arange(len(mapping))[matched],
            "fc1_row": mapping[matched],
        }
    )
    # The files a mapping belongs to, so a reader can refuse the wrong pair.
    sources = {"fc0_file": fc0_path.name, "fc1_file": fc1_path.name}
    table = table.assign(**sources)
    out_path = out_dir / f"event_matching_{tag}.parquet"
    table.to_parquet(out_path, index=False)

    calopart_path = out_dir / f"calopart_matching_{tag}.parquet"
    caloparts.assign(**sources).to_parquet(calopart_path, index=False)

    unmatched = pd.concat(
        [
            pd.DataFrame({"source_file": fc0_path.name, "row": unmatched_fc0}),
            pd.DataFrame({"source_file": fc1_path.name, "row": unmatched_fc1}),
        ],
        ignore_index=True,
    )
    unmatched_path = out_dir / f"unmatched_events_{tag}.parquet"
    unmatched.to_parquet(unmatched_path, index=False)

    print(f"wrote {out_path}  ({len(table)} rows)")
    print(f"wrote {calopart_path}  ({len(caloparts)} rows)")
    print(f"wrote {unmatched_path}  ({len(unmatched)} rows)")
    return table, caloparts.assign(**sources)


def read_mapping(name, tag, fc0, fc1):
    """`<name>_<tag>.parquet` from data/processed, refused unless it was built
    from this fc0/fc1 pair. `name` is `event_matching` or `calopart_matching`."""
    table = pd.read_parquet(output_paths.PROCESSED / f"{name}_{tag}.parquet")
    built = (table.fc0_file.iloc[0], table.fc1_file.iloc[0])
    if built != (Path(fc0).name, Path(fc1).name):
        raise SystemExit(f"{name}_{tag} was built from {built}, not ({fc0}, {fc1}); "
                         f"rerun: python scripts/analyze.py {fc0} {fc1} --steps match")
    return table
