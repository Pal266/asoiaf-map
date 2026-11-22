# Generating the resource layer (ore bodies and area resources)

## Inputs

- `lithology.npz` (rock, deposit, bedrock, province) from `world/lithology/scripts/lithology.py`
- `geology.npz` (age_ma, fault_type) from `world/geology/scripts/geology.py`
- `temperature.npz` and `precipitation.npz` (laterites)
- `soils.npz` (rsg) and `vegetation.npz` (biome): bog iron, peat, clay
- `canon_relief.npz` (feature_id) and `world/canon_relief/data/canon_features.json`: canon checks
- `data/resources.json`: deposit models (tract rules, grade-and-tonnage models and their sources), the density
  regression, area-resource rules, Earth reference counts

## Command

```
python3 world/resources/scripts/resources.py --lithology lithology.npz --geology geology.npz \
    --temperature temperature.npz --precipitation precipitation.npz --soils soils.npz \
    --vegetation vegetation.npz --canon-relief canon_relief.npz --out OUT
```

About 20 s; deterministic (seed in `resources.json`).

## Layer order

… → lithology → soils → vegetation → geology → **resources**.

## Steps

1. Permissive tracts of 21 deposit models from province, rock, age, nearness to intrusions, faults, evaporites or
   orogens, climate; placers on alluvium near lode-gold tracts. Outcrops closer than 30 km form one tract.
2. Number of deposits per tract from Singer's (2008) density regression (area and median tonnage); Poisson draw.
3. Tonnage from the lognormal grade-and-tonnage model (truncated at ±3σ); median grades; contained metal.
4. Area resources: coal measures, salt, bog iron, peat, clay, granite, limestone, marble, sandstone, slate,
   obsidian, sulfur.

## Outputs

- `resources.npz` [6, 1000, 1000]:
  - `deposit` int32: index into `deposits.json`, −1 where there is none; a cell with several deposits shows the
    first;
  - `permissive` uint32: one bit per deposit model, in the order of `resources.json`;
  - `area` uint16: one bit per area resource.
- `deposits.json`: models, area-resource names, and every deposit (type, cell, lat/lon, tonnage, grades, contained
  metal in tonnes)
- `stats.json`: counts and metal by type and region, permissive areas, area-resource areas, canon checks
- `preview.png`

## Expected values

- 69 525 deposits.
- Contained metal: copper 2.66 Gt, gold 395 000 t, silver 7.1 Mt, lead 7.9 Gt, zinc 18.5 Gt, tin 289 Mt,
  iron 1.24 × 10¹² t.
- Canon checks:
  - the Westerlands: 23 % of Westeros's gold and 56 % of its silver on 9.3 % of its area;
  - the Iron Islands: iron, lead, tin, and gold only 89 t;
  - the Fourteen Flames: gold, silver, copper, tin, iron.
