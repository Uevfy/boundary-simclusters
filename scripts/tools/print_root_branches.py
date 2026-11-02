"""Print branch names (and types) for the Events tree of one or more ROOT files.

Usage:
    python scripts/print_root_branches.py [file1.root file2.root ...]

With no arguments, scans all *.root files in the repo root.
"""

import sys
from collections import defaultdict
from pathlib import Path

import uproot

REPO_ROOT = Path(__file__).resolve().parents[2]
TREE_NAME = "Events"


def print_branches(path: Path) -> None:
    print(f"\n{'=' * 80}")
    print(f"{path.name}")
    print("=" * 80)

    with uproot.open(path) as f:
        if TREE_NAME not in f:
            print(f"  (no '{TREE_NAME}' tree found; top-level keys: {list(f.keys())})")
            return
        tree = f[TREE_NAME]

        by_prefix = defaultdict(list)
        for branch_name in tree.keys():
            prefix = branch_name.split("_", 1)[0] if "_" in branch_name else "(scalar)"
            by_prefix[prefix].append(branch_name)

        print(f"  {len(tree.keys())} branches, {tree.num_entries} entries\n")
        for prefix in sorted(by_prefix):
            names = by_prefix[prefix]
            print(f"  [{prefix}] ({len(names)})")
            for name in names:
                print(f"    {name}: {tree[name].typename}")


def main() -> None:
    args = sys.argv[1:]
    if args:
        paths = [Path(a) for a in args]
    else:
        paths = sorted(REPO_ROOT.glob("*.root"))

    if not paths:
        print("No ROOT files found.")
        return

    for path in paths:
        print_branches(path)


if __name__ == "__main__":
    main()
