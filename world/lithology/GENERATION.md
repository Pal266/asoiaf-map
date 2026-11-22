# Generating the lithology layer

## Inputs

- `history.npz` (province) from `world/history/scripts/history.py`
- `relief.npz` (elevation_m) from `world/relief/scripts/relief.py`
- `ice.npz` (sheet_frac) from `world/ice/scripts/ice.py`
- `hydrology.npz` (lake, discharge_m3s) from `world/hydrology/scripts/hydrology.py`
- `temperature.npz` and `precipitation.npz` from the climate layers run with `--ice ice.npz`
- `canon_relief.npz` (feature_id) from `world/canon_relief/scripts/canon_relief.py`: volcanic relief
- `data/lithology.json`: classes, province compositions, noise, deposit thresholds, GLiM shares

## Command

```
python3 world/lithology/scripts/lithology.py --history history.npz --relief relief.npz --ice ice.npz \
    --hydrology hydrology.npz --precipitation precipitation.npz --temperature temperature.npz \
    --canon-relief canon_relief.npz --out OUT
```

About 10 s; deterministic (fixed seed in `lithology.json`).

## Layer order

history → relief → ice ⇄ climate → hydrology → **lithology** → soils (`soils.py --lithology OUT/lithology.npz`),
geology (E16).

## Outputs

- `lithology.npz` [6, 1000, 1000]: `rock` (index into `order`), `deposit` (0 none, 1 alluvium, 2 evaporite,
  3 aeolian sand, 4 loess), `bedrock` (the rock under the deposits), `province` and `event_age_ma` (the
  history's provinces and event ages read through the warp, for the geology layer)
- `stats.json`: land shares against GLiM, deposits, shares by region
- `preview.png`: rock classes, equirectangular

## Expected values

Shares of land without water and ice: sediments 0.70 (GLiM 0.71), metamorphic 0.13 (0.14), plutonic 0.08 (0.08),
volcanic 0.08 (0.07), unconsolidated 0.22 (0.27); carbonates 26 % of sediments.
