# Generating the soil layer

## Inputs

- `relief.npz` (elevation_m) from `world/relief/scripts/relief.py`
- `ice.npz` (surface_m, sheet_frac) from `world/ice/scripts/ice.py`
- `temperature.npz` and `precipitation.npz` from the climate layers run with `--ice ice.npz`
- `hydrology.npz` (lake, discharge_m3s) from `world/hydrology/scripts/hydrology.py`
- `canon_relief.npz` (feature_id) from `world/canon_relief/scripts/canon_relief.py`, with
  `world/canon_relief/data/canon_features.json`: volcanic relief
- `data/soils.json`: the key's thresholds, the fertility of each group, Earth's relief distribution and the FAO areas
- for the Earth test only: the inputs of the ice and hydrology Earth tests, and `data/earth_soil_sites.json`

## Commands

```
python3 world/soils/scripts/earth_soil_validation.py --earth-data DIR                   # Earth test (optional)
python3 world/soils/scripts/soils.py --relief relief.npz --ice ice.npz --temperature temperature.npz \
    --precipitation precipitation.npz --hydrology hydrology.npz --canon-relief canon_relief.npz --out OUT
```

The world run takes about 10 s and the Earth test about 2 minutes. Both are deterministic.

## Steps (soils.py)

1. Thornthwaite potential evapotranspiration from the monthly temperature, then the moisture class and the number
   of dry months.
2. Local relief: height range over the 5 × 5 cells around each cell, mapped by rank onto Earth's distribution.
3. Volcanic parent material where the canon map draws volcanic relief; closed basins from the hydrology's terminal
   lakes, marshes and salt flats.
4. The key (`soilkey.classify`) gives each cell one Reference Soil Group; the fertility index follows from the group.

## Outputs

- `soils.npz` [6, 1000, 1000]: `rsg` uint8 (index into `groups` of `soils.json`), `fertility` uint8 (0–100)
- `stats.json`: areas per group, shares by region, the group at the canon places
- `preview.png`: soil groups, equirectangular

## Expected values

- **Earth test:** 52 % of 46 sites; area log RMSE 0.26 against the FAO estimates.
- **World:**
  - Westeros south of the Wall: Calcisols 19 %, Podzols 14 %, Luvisols 12 %; mean fertility 0.38;
  - beyond the Wall: Cryosols 94 %.
