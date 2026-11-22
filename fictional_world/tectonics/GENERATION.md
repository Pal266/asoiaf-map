# Tectonics — how to generate

All four steps are deterministic (seed in `data/tectonics.json`). Outputs go to a work directory and are not stored
(RULES.md rule 4).

```
W=/path/to/work
python3 scripts/supercontinent.py $W   # E01, ~4 s:  supercontinent.npz
python3 scripts/breakup.py $W          # E02, ~5 s:  breakup.npz
python3 scripts/drift.py $W            # E03, ~70 s: drift.npz, drift_snapshots.npz
python3 scripts/output.py $W           # E04, ~25 s: tectonics.npz (6 x 1000 x 1000), boundaries.npz, volcanoes.npz, tectonics_stats.json
```

Requires numpy, scipy and matplotlib. Each step writes a preview image to the work directory.
