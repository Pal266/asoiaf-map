# Generating the geology layer (rock ages, faults, seismic hazard)

## Inputs

- `plates.npz` (plate_cube) from `world/plates/scripts/plates.py`, with `world/plates/data/design.json` (boundary types)
- `history.npz` (province, event_age_ma, crust) from `world/history/scripts/history.py`, with
  `world/history/data/history.json` (spreading rates)
- `relief.npz` (elevation_m), `ice.npz` (sheet_frac)
- `lithology.npz` (rock, deposit, province, event_age_ma) from `world/lithology/scripts/lithology.py`
- `canon_relief.npz` (feature_id): the Fourteen Flames, for the Doom's age
- `data/geology.json`: periods, age rules, fault rules, hazard parameters
- for the Earth test only: PB2002 (`PB2002_steps.csv`, `PB2002_orogens.dig.txt` in the folder "original" of
  github.com/fraxen/tectonicplates) and `data/earth_hazard_sites.json`

## Commands

```
python3 world/geology/scripts/earth_hazard_validation.py --pb2002 DIR [--grid]     # Earth test (optional)
python3 world/geology/scripts/geology.py --plates plates.npz --history history.npz --relief relief.npz \
    --ice ice.npz --lithology lithology.npz --canon-relief canon_relief.npz --out OUT
```

The world run takes about 110 s (most of it the distance searches of the hazard) and the Earth test about 1 s. Both
are deterministic.

## Layer order

plates → history → relief → … → lithology → **geology** → resources.

## Steps (geology.py)

1. Rock ages: crystalline rocks from their province's event, sediments from the area-age power law within the
   limits of their setting, volcanics as young as their arc or rift, deposits Quaternary; periods of the ICS chart.
2. Faults: plate contacts (subduction, ridges, transforms) and province contacts on continental crust (thrusts,
   normal faults, sutures), each with the age of its last movement and an activity class.
3. Seismic hazard (hazard.py): plate boundaries by class and velocity, active and dormant deforming zones; GSHAP
   classes.

## Outputs

- `geology.npz` [6, 1000, 1000]:

  | Array | Type | Content |
  |---|---|---|
  | `age_ma` | float32 | formation age of the surface rock, Ma (NaN on water and ice) |
  | `period` | uint8 | index into `periods.names` (255 none) |
  | `fault_type` | uint8 | 0 none, 1 subduction megathrust, 2 spreading ridge, 3 transform, 4 thrust, 5 normal, 6 suture |
  | `fault_age_ma` | float32 | age of the fault's last movement (NaN where no fault) |
  | `fault_activity` | uint8 | 0 none, 1 active, 2 dormant, 3 inactive |
  | `hazard` | uint8 | hazard index × 100 |
  | `hazard_class` | uint8 | 0 low, 1 moderate, 2 high, 3 very high (GSHAP) |

- `stats.json`: period and era shares, area-age slope, fault lengths, hazard by region, the canon places
- `preview_age.png` (periods, faults in black), `preview_hazard.png`

## Expected values

- **Earth test:** 87 % of 71 cities in their GSHAP class (leave-one-out 85 %).
- **World:**
  - eras of the land: Cenozoic 56 %, Mesozoic 9 %, Paleozoic 12 %, Precambrian 23 %; area-age slope −0.97;
  - hazard on land: low 70 %, moderate 9 %, high 19 %, very high 2 %; Westeros 72 / 18 / 4 / 5 %.
