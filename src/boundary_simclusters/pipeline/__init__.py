"""The standard run, one module per step, each with `run(sample, ...)`:

    match_events -> paired_calopart_summary -> make_figures -> classify_surfaces
        -> classify_decompositions

Run from the command line by scripts/analyze.py (`--steps` for a subset), and
all at once by `numbers.measure`.
"""
