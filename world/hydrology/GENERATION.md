# Generating the hydrology layer

## Inputs

- `relief.npz` (elevation_m) from `world/relief/scripts/relief.py`
- `ice.npz` (surface_m, sheet_frac) from `world/ice/scripts/ice.py`
- `temperature.npz` from `world/climate/scripts/climate.py --ice ice.npz`
- `precipitation.npz` from `world/precipitation/scripts/precipitation.py --ice ice.npz`
- `canon_relief.npz` (feature_id) from `world/canon_relief/scripts/canon_relief.py`, for the checks
- map M1 (`world/canon_outline/research/SOURCES.md`), for the canon rivers and lakes
- `data/hydrology.json` (parameters) and `data/canon_rivers.json` (canon rivers and lakes traced on M1)
- for the Earth test only: ETOPO1 10′, Natural Earth glaciated areas and 1:50m rivers (downloaded once into
  `--earth-data`), and `data/earth_rivers.json`

## Commands

```
python3 world/hydrology/scripts/earth_hydro_validation.py --earth-data DIR             # Earth test (optional)
python3 world/hydrology/scripts/canon_rivers.py --map M1.jpg --out OUT
python3 world/hydrology/scripts/hydrology.py --relief relief.npz --ice ice.npz --temperature temperature.npz \
    --precipitation precipitation.npz --canon-relief canon_relief.npz --canon-rivers OUT/canon_rivers.npz --out OUT
```

`canon_rivers.py` takes about 2 minutes, `hydrology.py` about 30 s and the Earth test about 2 minutes. All are
deterministic.

## Steps

1. **Canon rivers and lakes** (`canon_rivers.py`):
   - each traced river is snapped to the line drawn on M1 (least-cost path through dark pixels in an 18 px
     corridor);
   - the water inside each lake envelope is taken as lake;
   - both are placed on the 10 km cells.
2. **Water balance:** Thornthwaite potential evapotranspiration × 1.5, then runoff from the Budyko curve
   (Fu, w = 2.6). On ice sheets all precipitation runs off.
3. **Routing surface:** the ice surface (bedrock where there is no ice), with the canon rivers carved 500 m deep,
   falling from source to mouth, and continued to the coast.
4. **Routing:** priority-flood filling from the sea, steepest descent over 8 neighbours, accumulation downstream
   (never below zero).
5. **Lakes** on the uncarved surface:
   - depressions ≥ 1000 km² and ≥ 10 m deep at their deepest;
   - overflowing ones keep water deeper than 50 m, and their shallower parts are marsh;
   - terminal ones shrink until their open-water evaporation equals the inflow; up to 20 m above the lake their
     floor is salt flat;
   - the canon lakes are lakes, terminal when nothing flows out.
6. **Strahler order** of channels with ≥ 10 m³/s.

## Outputs

- `hydrology.npz`, all [6, 1000, 1000]:

  | Field | Type | Content |
  |---|---|---|
  | `discharge_m3s` | float32 | river discharge |
  | `runoff_mm` | float32 | annual runoff |
  | `downstream` | int8 | receiving neighbour 0–7; −1 at sea and at terminal sinks |
  | `lake` | uint8 | 1 lake, 2 terminal lake, 3 marsh, 4 salt flat |
  | `strahler` | uint8 | stream order |

- `stats.json`: runoff, lake areas, largest rivers, canon rivers and lakes.
- `preview.png`.
- From `canon_rivers.py`: `canon_rivers.npz` (river and lake cells, river cell sequences), `canon_rivers.json`,
  `check.png`.

## Expected values

- **Earth test:**
  - runoff to the sea 40 300 km³/yr;
  - river discharge log RMSE 0.82, limited by the precipitation model;
  - drainage areas log RMSE 0.36.
- **World:**
  - land runoff 298 mm/yr; 38 300 km³/yr to the sea;
  - Rhoyne at Volantis 1 160 m³/s, Trident 3 120, White Knife 3 990, Mander 240;
  - Bleeding Sea and Shrinking Sea terminal.
