# Generating the temperature layer

## Inputs

- `relief.npz` (elevation_m) from `world/relief/scripts/relief.py`
- `world/climate/data/climate.json`: model parameters calibrated on Earth, the station correction and the world's
  warm-drift box
- for the Earth test only: `world/climate/data/earth_stations.json` and the pip package `global-land-mask`

## Commands

```
python3 world/climate/scripts/earth_validation.py --out EARTH       # test and refit on Earth (optional)
python3 world/climate/scripts/climate.py --relief relief.npz --out OUT
```

The world run takes about 30 s and the Earth test about 30 s. Both are deterministic.

## Steps (climate.py)

1. **Land fraction** of every 2° cell, from the 10 km land mask (elevation > 0), sampled at 4 × 4 points per cell.
2. **Two-layer energy balance model** (`ebm.run_two_layer`), 12 years of spin-up, daily steps:
   - insolation from Earth's orbit;
   - planetary albedo with snow and ice;
   - surface and air layers, with stability-dependent exchange and window emission;
   - diffusion and zonal-wind advection of the air;
   - continentality from the distance to the sea upwind;
   - zero-layer sea ice;
   - the warm drift (Q-flux) in the box `world_qflux_box`.

   Output: monthly mean sea-level surface temperature on the 2° grid.
3. **To the 10 km grid:** bilinear interpolation, minus 6.5 K/km × ground elevation, minus the statistical
   correction on land: b0 + b1·min(z, 2.5 km) + b2·M. Its coefficients move from winter to summer with the month and
   the hemisphere and fade over 15° at the equator.
4. **Statistics:** places, the Haunted Forest's warmest month, the Land of Always Winter.

## Outputs

- `temperature.npz`: `t_month_dC` int16 [12, 6, 1000, 1000] (tenths of °C), `sea_ice` uint8 [12, 90, 180]
- `stats.json`, `preview.png` (July over January)

## Expected values

- **Earth test:** RMSE 5.4 / 4.0 / 3.6 °C for January / July / year at 78 stations; after correction,
  leave-one-out 5.0 (winter) and 3.2 °C (summer).
- **World:**
  - global mean 15.6 °C;
  - Winterfell −9.0 / 11.3 °C in January / July;
  - Haunted Forest warmest month 7.6 °C.
