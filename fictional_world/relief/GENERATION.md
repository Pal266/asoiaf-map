# Relief — how to generate

Requires the tectonics layer in the same work directory (`../tectonics/GENERATION.md`).

```
W=/path/to/work
python3 scripts/relief.py $W      # E01, ~2 min: relief.npz (elevation_m, subgrid_relief_m, sediment_m), relief_stats.json, E01_relief.png, E01_relief_zoom.png
```

Deterministic (seed in `data/relief.json`). Outputs are not stored (RULES.md rule 4).
