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
2. Continental crust = land, or ocean within the margin width of land (300 km; 60 km where a trench is closer than 300 km).
3. Ocean floor age: among the ridges bounding the cell's plate (all plates for a ridge listed with plate 0) and within their reach (half-rate × opening age), the nearest; age = distance / half-rate. Extinct ridges (`fossil_ridges`) inside a plate give end age + distance / half-rate, up to their start age. Floor no ridge reaches takes 160 Ma.
4. Oceanic plateaus: a seeded fractal field (seed + 1) on ocean cells with floor older than 30 Ma and more than 300 km from trenches and land; the highest values covering 5.11 % of the ocean area become province 9.
5. Basement depth from age, GDH1: 2600 + 365√t (t < 20 Myr), 5651 − 2473·exp(−0.0278 t) otherwise.
6. Provinces on land, in order (later steps override earlier ones):
   - blocks: polygon test on warped coordinates (fractal, 4 octaves 4000–500 km, Hurst 0.7, 3°), same plate only; blocks with `belts` split at the 45% quantile of a second fractal field into craton cores and 600 Ma belts;
   - land left over: platform, 600 Ma;
   - margins: cells of the given plate within the given distance of its coast that faces the boundary (coast cells within 1000 km of the boundary line);
   - canon features: province and age from E05.
7. Shelf cells copy the province and age of the nearest land cell; ocean cells get province 0 and their floor age.

## Outputs

- `history.npz`: crust, province (codes 0–9 in history.json), event_age_ma, ocean_depth_m — all [6, 1000, 1000]
- `stats.json`: province areas and checks (including the floor share by spreading class) (continental share, floor age, depth, narrow sea age, canon features in the wrong province)
- `preview.png`: provinces and ocean floor age

## Checks and expected values

continental crust 40.1%; oceanic plateaus 5.11 % of the ocean; floor age max 160 Ma, mean 60 Ma; mean basement depth 4834 m; floor by spreading class slow 20 %, intermediate 27 %, fast 47 %, older than every ridge 6 %; narrow sea floor ≤ 25 Ma (23.4); narrow sea floor ≤ 25 Ma; no canon feature in the wrong province.
