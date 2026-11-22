# Generating the vegetation layer

## Inputs

- `relief.npz` (elevation_m) from `world/relief/scripts/relief.py`
- `ice.npz` (surface_m, sheet_frac) from `world/ice/scripts/ice.py`
- `temperature.npz` and `precipitation.npz` from the climate layers run with `--ice ice.npz`
- `hydrology.npz` (lake, discharge_m3s) from `world/hydrology/scripts/hydrology.py`
- `soils.npz` (rsg) from `world/soils/scripts/soils.py`, run with `--lithology`: peat soils (Histosols)
- `canon_relief.npz` (feature_id) from `world/canon_relief/scripts/canon_relief.py`: canon checks only
- `data/vegetation.json`: biomes, the key's limits, colours, the mapping to the WWF biomes and their Earth areas
- for the Earth test only: the Earth data of the ice and hydrology tests, and `data/earth_veg_sites.json`

## Commands

```
python3 world/vegetation/scripts/earth_veg_validation.py --earth-data DIR                  # Earth test (optional)
python3 world/vegetation/scripts/vegetation.py --relief relief.npz --ice ice.npz --temperature temperature.npz \
    --precipitation precipitation.npz --hydrology hydrology.npz --soils soils.npz \
    --canon-relief canon_relief.npz --out OUT
```

The world run takes about 15 s and the Earth test about 70 s. Both are deterministic.

## Layer order

history → relief → ice ⇄ climate → hydrology → lithology → soils → **vegetation**.

## Steps (vegetation.py, through vegkey.py)

1. Thornthwaite potential evapotranspiration; a 150 mm bucket gives the actual evapotranspiration and the moisture
   index α; growing degree-days, coldest and warmest month, Köppen's Cs test, dry months.
2. The key: polar desert, tundra and alpine by summer heat; desert and xeric shrubland by α; the tropics by winter
   cold (savanna, dry forest, rainforest); Mediterranean woodland; the boreal zone (forest-tundra, taiga);
   grassland; subtropical, conifer, mixed and broadleaf forests.
3. Open land above 2500 m (below 55° latitude) is alpine.
4. Marshes and peat soils are wetlands (swamp forest in forest climates); large rivers through dry land carry
   riparian woodland; salt flats are bare.
5. Net primary productivity from the Miami model.

## Outputs

- `vegetation.npz` [6, 1000, 1000]: `biome` uint8 (index into `biomes` of `vegetation.json`), `npp_gm2` uint16
  (g dry matter / m² / yr)
- `stats.json`: areas, shares by region, Earth-equivalent (WWF) shares against Earth's, canon checks (Haunted
  Forest, Dothraki sea, canon deserts and the Neck), the biome at the canon places, global NPP
- `preview.png`: biomes, equirectangular

## Expected values

- **Earth test:** 64 % of 78 sites; area log RMSE 0.14 against the WWF biomes; NPP 115 Pg/yr.
- **World:**
  - Westeros south of the Wall: Mediterranean woodland 16 %, temperate broadleaf forest 16 %, xeric shrubland 11 %;
  - beyond the Wall: tundra 90 %;
  - Dothraki sea 60 % grassland; Red Waste 75 % desert;
  - global NPP 103 Pg/yr.
