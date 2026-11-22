# Generating the precipitation layer

## Inputs

- `relief.npz` (elevation_m) from `world/relief/scripts/relief.py`
- `canon_relief.npz` (feature_id) from `world/canon_relief/scripts/canon_relief.py`, for the desert checks
- `world/climate/data/climate.json`: the temperature model is rerun inside (same run as the temperature layer)
- `world/precipitation/data/precipitation.json`: model parameters calibrated on Earth
- for the Earth test only: `data/earth_precip_stations.json`, `data/earth_zonal_gpcp.json`, the pip package
  `global-land-mask`, and ETOPO1 10′ topography (downloaded once into `--earth-data`)

## Commands

```
python3 world/precipitation/scripts/earth_precip_validation.py --earth-data DIR --out EARTH   # Earth test (optional)
python3 world/precipitation/scripts/precipitation.py --relief relief.npz --canon-relief canon_relief.npz --out OUT
```

The world run takes about 1 minute. Both runs are deterministic.

## Steps (precipitation.py)

1. **Temperature and sea ice:** the two-layer energy balance model of the temperature layer, monthly, sea level,
   2°.
2. **Land fraction and mean ground height** of every 2° cell, from the 10 km relief.
3. **Moisture model** (`moisture.run`), two passes of the annual cycle; each month:
   - winds: zonal belts, meridional cells about the ITCZ, thermal and geostrophic components, the Gill monsoon
     circulation;
   - evaporation from sea (wind, temperature) and land (limited by last pass's rain);
   - precipitation rate from convergence, upslope wind, Hadley subsidence and Sverdrup vertical motion;
   - one sparse solve of div(WV) − div(K∇W) + kW = E; P = kW.
4. **To the 10 km grid:** bilinear interpolation of the monthly totals.
5. **Statistics:** places, means over the canon deserts, Westeros' arid share.

The model must run at 2°, the resolution it was calibrated at (upwind numerical diffusion depends on the grid step).

## Outputs

- `precipitation.npz`: `p_month_mm` uint16 [12, 6, 1000, 1000] (mm per month)
- `stats.json`, `preview.png` (annual precipitation)

## Expected values

- **Earth test:** log10 RMSE 0.39 at 109 stations, 56 % within a factor 2; arid-station bias +0.22 log; global
  mean 1026 mm.
- **World:** global mean 985 mm; Winterfell 690 mm; Dornish desert median 268 mm; Red Waste median 85 mm.
