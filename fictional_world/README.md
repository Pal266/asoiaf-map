# Fictional world

A second world, built from scratch: no canon map, no canon text. It follows the same rules as the rest of the
repository (`../RULES.md`). The owner's design decisions, which take the place of canon here, are recorded in
`DECISIONS.md`.

## Layers

Each level is built from the levels above it. Every layer has its own folder with research, experiments,
generation scripts and documents (Rules 4 and 5); generated layers are not stored.

| # | Layer | Folder | Status |
|---|---|---|---|
| 1 | Planet: Earth's parameters | `planet/` | done |
| 2 | Global grid: cubed sphere, 6 × 1000 × 1000 cells, ~10 km | `grid/` | done (taken from the first world) |
| 3 | Tectonic history: plates, the breakup of a supercontinent, crust, ocean floor ages, tectonic provinces | `tectonics/` | done (E01–E05) |
| 4 | Relief | `relief/` | done (E01) |
| 5 | Ice and temperature | `ice/`, `climate/` | |
| 6 | Precipitation | `precipitation/` | |
| 7 | Hydrology | `hydrology/` | |
| 8 | Lithology | `lithology/` | |
| 9 | Soils | `soils/` | |
| 10 | Vegetation | `vegetation/` | |
| 11 | Rock ages, faults, seismic hazard | `geology/` | |
| 12 | Mineral resources | `resources/` | |
| 13 | Game region | `region/` | chosen after the climate (decision 4) |

Layers 4–12 can reuse the models of the first world (`world/` on branch `world-geological-history`), which were
tested on Earth. What changes is layer 3: there the first world started from the canon map, while this world starts
from a simulated plate history.
