# Generating the ice layer

## Inputs

- `relief.npz` (elevation_m, subgrid_relief_m) from `world/relief/scripts/relief.py`
- `world/ice/data/ice.json`: mass balance, glacier and ice-sheet parameters, the start state
- `world/climate/data/climate.json`, `world/precipitation/data/precipitation.json`: the climate models are run inside
- for the Earth test only: ETOPO1 10′ and Natural Earth 1:10m glaciated areas, downloaded once into `--earth-data`

## Commands

```
python3 world/ice/scripts/earth_ice_validation.py --earth-data DIR --out EARTH      # Earth test (optional)
python3 world/ice/scripts/ice.py --relief relief.npz --out OUT
```

The world run takes about 25 minutes, mostly in the climate models, which run once per pass. It is deterministic.

## Order of the climate layers

The ice and the climate depend on each other, so the ice layer finds them together from the relief. The temperature
and precipitation layers are then generated with the ice:

```
relief → ice → temperature (climate.py --ice ice.npz) → precipitation (precipitation.py --ice ice.npz)
```

## Steps (ice.py)

1. **Start:** the south polar continent with ice sheets wherever its bare-ground annual mean is below −5 °C (its
   geological history keeps it at the pole); everything else bare.
2. **Climate of the current surface** (`world_climate.py`): the temperature model sees the land height (snow line)
   and the ice sheets (albedo); the moisture model follows.
3. **Inner loop**, with that climate:
   - surface mass balance of every cell at the ice surface (degree days; colder and drier on raised ice);
   - ice sheets: positive-balance areas ≥ 5000 km², grown into the ablation zone while their balance minus calving
     stays positive;
   - plastic thickness inwards from the margin (τ0 75 kPa), with local isostasy.
4. **Glaciers** outside the sheets, from each cell's sub-grid relief.
5. Repeat from step 2 until the sheets change by less than 50 000 km².

## Outputs

- `ice.npz`, all [6, 1000, 1000]:
  - `ice_thickness_m` float32;
  - `surface_m` float32: ice surface on land, depressed bedrock under the ice, sea floor at sea;
  - `sheet_frac` uint8: 255 = ice sheet, otherwise glacier share × 254.
- `stats.json`: passes, areas, volume, the Land of Always Winter.
- `preview.png`: ice over shaded relief.

## Expected values

- **Earth test:** today's surface, CSI 0.64 (σ 5.5); coupled from bare ground (σ 3.5), Antarctica regrows,
  Greenland does not.
- **World:**
  - ice sheets 14.74 M km², 13.99 of them on the south polar continent, max 6.6 km;
  - Land of Always Winter 2.4 % ice sheet, a 28 000 km² mountain ice cap;
  - no ice in Westeros south of 70° N.
