# E06. Present-day plates and the land frame of the planet

## Goal

1. Draw the present-day plates of `research/E05_tectonic_interpretation.md` on the globe.
2. Place the lands outside the map as frames with the target areas of `world/canon/CONSTRAINTS.md` (decision 6,
   revised), so that land covers ~29 % of the planet (decision 4).
3. Check that the plates explain the canon features they were designed for.

## Input data

- The canon outline layer (`world/canon_outline/GENERATION.md`).
- `data/design.json`: unknown-land polygons and target areas, coast parameters, plate seeds, boundary lines.
- Script `scripts/plates.py` (about 5 s).

## Method

- **Design grid:** equirectangular, 0.5° (720 × 360), the coarse grid that `world/grid/GRID.md` allows for tectonics;
  results are also sampled onto the global cubed sphere.
- **Unknown lands:** each frame is a polygon (or, for the south polar continent, a cap). Its coast is the level set of
  (signed distance to the polygon + a fractal field), the field being a sum of plane waves along random directions on
  the sphere (wavelengths 3 000–375 km, Hurst exponent 0.7, amplitude 450 km, fixed seed). The level is found by
  bisection so that the frame's land matches its target area. Within 250 km of the map's edge, unknown cells continue
  the nearest canon coast (land where the map's edge is land, sea where it is sea), so coasts run on across the edge.
- **Plates:** boundaries are polylines of great-circle segments (the polar ring runs along 62° S in 5° steps);
  plates are the regions between them, flood-filled from one seed per plate, with wrap-around in longitude.

## Iterations

| Run | Change | Result |
|---|---|---|
| a | First design | Areas below target: land continued across the map's edge was counted for the frame but stored as canon |
| b | Continued coasts counted with their frame | All targets met; Sothoryos and Ulthos overlapped across the Sothoryos–Ulthos ridge (a ridge may not cut a continent) |
| c | Ulthos frame moved east, ridge to 70.5° E, trench nearer Valyria | The south polar land crossed the ring ridge: great-circle segments between distant vertices bulge towards the pole |
| d | Polar ring drawn along 62° S in 5° steps; polar land kept south of 64° S | Final run |

## Results (run d)

**Land:** 147.1 million km², **28.8 %** of the planet (Earth 29 %). Canon land on the design grid: 32.1 million km².

| Frame | Target | Result |
|---|---|---|
| Essos east | 20.0 M km² | 20.0 |
| Sothoryos south | 31.0 | 31.0 |
| Ulthos | 10.0 | 10.0 |
| Southern continent | 39.0 | 39.0 |
| South polar continent | 14.0 | 14.0 |
| Land of Always Winter north of the map | 1.0 | 1.0 |

**Plates:**

| Plate | Area, M km² | Land, M km² |
|---|---|---|
| Westeros | 65.1 | 8.1 |
| Essos | 83.2 | 42.0 |
| Sothoryos | 80.5 | 32.9 |
| Ulthos | 47.5 | 10.8 |
| Sunset (oceanic) | 126.5 | 0.4 |
| Southern continent | 78.4 | 38.9 |
| South polar | 28.8 | 14.0 |

**Checks:**

- Places lie on their designed plates: Winterfell, King's Landing, Sunspear and Dragonmont on Westeros; Volantis,
  Valyria, Qarth, Asshai, Ibben and Great Moraq on Essos; the Summer Isles and the Basilisk Isles on Sothoryos.
- Valyria lies 244 km from the Summer Sea trench, an arc–trench distance like the Hellenic arc (Santorini lies
  ~120 km north of Crete, which itself lies north of the trench).
- Dragonmont lies 358 km from the narrow-sea ridge, on the rift flank (volcanism away from the axis of a young rift).
- No land mass is cut by a ridge or a transform except for slivers under 0.1 % of their area; land on two plates
  occurs only along trenches and collision zones: Essos 0.7 % beyond the East Essos trench, Ulthos 2 % across the Jade
  arc collision, as forearc and accreted slivers.
- Two runs give byte-identical outputs.

## Conclusions

- **Goals 1–3 reached.** The plates of E05 are drawn, the unknown lands have frames with their target areas, and the
  plates contain the canon features they were designed to explain.
- The frames are not final coasts: the geological history will shape the unknown lands (stage 3); the frames give
  their places and sizes. The boundaries are schematic polylines, to be refined by the history.
- The Land of Always Winter frame reaches near the north pole, forming land under the polar cap north of the
  Shivering Sea, which canon describes as a polar cap.
