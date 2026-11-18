# Generating the geological history layer

## Inputs

- `plates.npz` from `world/plates/scripts/plates.py` (plate_cube, land_cube)
- `canon_relief.npz` from `world/canon_relief/scripts/canon_relief.py` (feature_id)
- `world/plates/data/design.json` (boundaries) and `world/history/data/history.json` (scenario, ridges, blocks, margins, canon features, seed 267)

## Command

```
python3 world/history/scripts/history.py --plates plates.npz --relief canon_relief.npz --out OUT
```

About 4 minutes on one CPU. The same inputs give the same output.

## Steps

1. Cell centres of the 6×1000×1000 cubed sphere; distances to every boundary line (KD-tree on unit vectors, cached).
2. Continental crust = land, or ocean within the margin width of land (250 km; 60 km where a trench is closer than 300 km).
3. Ocean floor age: the nearest ridge bounding the cell's plate; age = distance / half-rate, capped at the ridge's opening age. Default 120 Ma if a plate has no ridge.
4. Basement depth from age, GDH1: 2600 + 365√t (t < 20 Myr), 5651 − 2473·exp(−0.0278 t) otherwise.
5. Provinces on land, in order (later steps override earlier ones):
   - blocks: polygon test on warped coordinates (fractal, 4 octaves 4000–500 km, Hurst 0.7, 3°), same plate only; blocks with `belts` split at the 45% quantile of a second fractal field into craton cores and 600 Ma belts;
   - land left over: platform, 600 Ma;
   - margins: cells of the given plate within the given distance of a boundary;
   - canon features: province and age from E05.
6. Shelf cells copy the province and age of the nearest land cell; ocean cells get province 0 and their floor age.

## Outputs

- `history.npz`: crust, province, event_age_ma, ocean_depth_m — all [6, 1000, 1000]
- `stats.json`: province areas and checks (continental share, floor age, depth, narrow sea age, canon features in the wrong province)
- `preview.png`: provinces and ocean floor age

## Checks and expected values

continental crust 38.7%; floor age max 180 Ma, mean 80 Ma; mean basement depth 5023 m; narrow sea floor ≤ 25 Ma; no canon feature in the wrong province.
