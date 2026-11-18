# Layer "plates": how it is generated

The present-day tectonic plates and the land frame of the whole planet (canon land plus frames of the unknown lands).
The layer is not stored in the repository (Rule 4).

## What the layer is

| File | Content |
|---|---|
| `plates.npz` → `land`, `source`, `plate` | design grid, equirectangular 0.5° (720 × 360, rows north to south, columns from 180° W): land (1/0); source of the land (0 sea, 1 canon map, 2 + k unknown frame k of `design.json`); plate id |
| `plates.npz` → `plate_cube`, `land_cube` | the same plate id and land on the global cubed sphere `[6, 1000, 1000]` (nearest design cell) |
| `plates.geojson` | boundary lines with name, type (ridge, trench, transform) and the plates they separate |
| `stats.json` | areas per frame and plate, checks |
| `preview.png` | plates, land (canon land lighter) and boundaries (ridge yellow, trench red, transform cyan), 1 440 × 720 |

## Input

- `canon_outline.npz` from `world/canon_outline/GENERATION.md`, step 2.
- `data/design.json` (in the repository): frames of the unknown lands (polygons or a polar cap, target areas), coast
  parameters, plate seeds and boundary lines in longitude/latitude.

## Steps

From `world/plates/scripts/`:

```
python3 plates.py --outline <work>/canon_outline/canon_outline.npz --out <work>/plates
```

About 5 s.

1. Canon land: each design cell takes `land_frac` and `known_frac` of the canon outline at its centre; cells with
   `known_frac` > 0.5 are canon, land where `land_frac` > 0.5.
2. Next to the map's edge (250 km), unknown cells take the land/sea state of the nearest canon edge cell.
3. Each frame: signed distance to its polygon (positive inside; for the cap, distance south of its latitude) plus a
   fractal field (24 plane waves per octave along random directions, wavelengths 3 000, 1 500, 750, 375 km, Hurst
   exponent 0.7, amplitude 450 km, seed 266). The threshold is found by bisection so that the frame's land,
   including its share of step 2, equals the target area.
4. Plates: boundary lines are densified along great circles (0.1°) and rasterised; the regions between them are
   labelled with wrap-around in longitude and named by the plate seeds; boundary cells join the nearest plate.
5. Checks: total land share, frame areas, the plate of named places, Valyria's distance from the Summer Sea trench,
   Dragonmont's distance from the narrow-sea ridge, land masses lying on more than one plate.

Why: `research/E05_tectonic_interpretation.md` and `experiments/E06_plate_design.md`.

## Checks (run d of E06)

| Check | Expected | Result |
|---|---|---|
| Land share | ~29 % | 28.8 % |
| Frame areas | targets | all within 0.1 % |
| Named places | on their plates | all |
| Land cut by ridges or transforms | none | slivers < 0.1 % |
| Determinism | identical | identical |
