# Generating the relief layer

## Inputs

- `plates.npz` (plates layer), `canon_relief.npz` (canon relief layer), `history.npz` (history layer)
- `world/plates/data/design.json`, `world/relief/data/relief.json` (all parameters; seed 268)

## Command

```
python3 world/relief/scripts/relief.py --plates plates.npz --relief canon_relief.npz --history history.npz --out OUT
```

About 5 minutes on one CPU. Deterministic.

## Steps

1. Cell centres; seam-aware neighbours (`cubesphere.halo_index`). Four fractal fields from the seed:
   - plain noise;
   - range noise;
   - sea-floor noise;
   - ridged noise 1 − |f|.
2. Distances: to trenches, from land to sea, from sea to land (KD-tree on unit vectors). A cell is "active" if a
   trench is within 300 km.
3. Plain: province level and amplitude (`plain_m`), each blurred by 40 box passes (σ ≈ 50 km). The value is
   level + amplitude × noise, at least ¼ of the level. Toward the coast it is ramped by (d/W)^0.7, where
   W = 150 km (passive) or 40 km (active).
4. Ranges, each giving a crest height C and a profile s (0 at the edge, 1 on the crest):
   - generic ranges by province;
   - active-margin ranges: the plate's land within `width_km` of the coast facing the boundary;
   - canon ranges: C × (0.7 + 0.3 × drawn relief density).

   Then:
   - range term = C × sin(π s / 2) × (0.8 + 0.25 × noise) × (0.5 + 0.8 × ridged);
   - blurred by 10 passes;
   - zero at the coastline cell;
   - maximum with the cones of point features (Marahai: its ring of land within 200 km takes a range profile).
5. Canyon terms (sine profile, blurred) are subtracted. The Neck is pulled toward 30 m with a blurred interior
   weight. Land is kept at 1 m or higher.
6. Sea floor: − GDH1 depth + 0.57 × sediment (3.25 m/Myr × age up to 1 km, plus 3134 m × exp(−(d/554 km)²) from
   land; calibrated to GlobSed) − trench term + abyssal noise + 2000 m × blurred oceanic-plateau
   mask (province 9) + residual depth (fractal, rms 261 m, 4000–500 km, set so that 6.3 % of the ocean is deeper than 5.5 km). On continental-margin cells the
   floor is a shelf to 130 m over the inner 40 % of the margin, then a smoothstep slope to the nearest oceanic
   floor (margin width from history.json). Stepstones are set to −40 m. Sea is kept at −1 m or lower.
7. Sub-grid relief = 40 + 0.9 × range term + 0.5 × canyon term on land; 300 m at sea.

## Outputs

- `relief.npz`: `elevation_m`, `subgrid_relief_m`, `sediment_m` (sea floor) — float32 [6, 1000, 1000]
- `stats.json`: sediment budget against GlobSed, hypsometry, coast consistency, highest cell of every canon range against its analogue
- `preview.png`: shaded elevation, equirectangular 1440 × 720

## Expected values (run p, with E07 run g)

- mean land 723 m; land above 1000 m 19 %, above 2000 m 4.1 %;
- mean ocean depth 3880 m, median 4365 m; deeper than 5.5 km 5.6 %;
- sediment: margin zone 3044 m, far ocean 404 m, total 337.1 million km³;
- highest cell 5823 m, deepest −8795 m;
- `coast_consistent` true.
