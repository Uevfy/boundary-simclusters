"""The standard run, one module per step, each with `run(sample, ...)`:

    match_events -> paired_calopart_summary -> make_figures -> classify_surfaces
        -> classify_decompositions

Run from the command line by scripts/analyze.py (the first four, nanoML only;
`--steps` for a subset) and scripts/classify.py (classify_decompositions, which
needs RECO as well), and all at once by `numbers.measure`.
"""
