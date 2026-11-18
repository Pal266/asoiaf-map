# Layer "canon relief": how it is generated

The relief features drawn on the official map: envelopes of mountain ranges, hills, volcanic areas, deserts, basins,
inland lakes and wetlands, as targets for the geological history. The layer is not stored in the repository
(Rule 4); these steps generate it.

## What the layer is

| File | Content |
|---|---|
| `canon_relief.npz` → `feature_id` | int16 `[6, 1000, 1000]`: id of the canon feature at the cell centre (0 = none); ids as in `data/canon_features.json` |
| `canon_relief.npz` → `relief` | float32 `[6, 1000, 1000]`: density of drawn relief strokes inside mountain, hill and volcanic envelopes, 0–1, mean of 4 × 4 samples per cell; 0 elsewhere |
| `canon_features.geojson` | the envelopes (polygons) and points in longitude/latitude, with name, type, region, canon fact, source, area and latitudes |
| `stats.json` | measures per feature and per type |
| `preview_map.png` | envelopes (numbered) and strong relief over M1, 2 262 × 1 500 px |

Grid and placement as in `world/canon_outline/GENERATION.md`.

## Input

- Map M1 (`world/canon_outline/research/SOURCES.md`).
- The land mask: `land_mask.png` from step 1 of `world/canon_outline/GENERATION.md`.
- `data/canon_features.json` (in the repository): 34 features drawn by eye on gridded views of M1, in map pixels,
  each with type, region, canon fact and source. This file is research data: to change a feature, edit it there.

## Steps

From `world/canon_relief/scripts/`:

```
python3 canon_relief.py --map <path to M1> --mask <work>/mask/land_mask.png --out <work>/canon_relief
```

About 45 s.

1. Each envelope is filled on the map; for mountains, hills, deserts, volcanic areas, wetlands and basins it is
   clipped to land; lakes and the sunken land bridge are not.
2. Relief strokes: pixels at least 0.12 darker, or 0.10 brighter and weakly saturated (< 0.25), than the local
   background (Gaussian, σ = 15 px), excluding rivers, roads, labels and icons (dilated 3 px), a 10 px coastal band,
   the fold lines (rows 1985 and 3986 ± 20 px, columns 1650 and 7240 ± 12 px) and forests. Their density is the
   share of stroke pixels among usable land pixels (Gaussian, σ = 15 px), where at least half the neighbourhood is
   usable.
3. Each grid cell takes the feature at its centre and the mean stroke density of 4 × 4 sample points, using the
   placement of the canon outline (map pixel lookup at λ, λ + 360°, λ − 360°).

Why: experiment `experiments/E04_canon_relief.md`.

## Checks

| Check | Result |
|---|---|
| Area of each feature, grid vs map (equal-area placement) | within 3 % for every feature, most within 0.5 % |
| Determinism | two runs give byte-identical `canon_relief.npz` and `canon_features.geojson` |
| Visual check | envelopes and strong relief over M1 in `preview_map.png` |

## Known limits

- Envelopes are generous by design and drawn by eye; their edges are not canon.
- The stroke density is weak where relief is drawn in brown on orange land (Red Mountains) or with light strokes; it
  measures drawing, not height.
