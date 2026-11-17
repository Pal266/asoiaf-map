# E01. Land/sea mask from the official map

## Goal

1. Separate land from sea on map M1 (`research/SOURCES.md`) automatically, so that the canon outlines can be used.
2. Find the map's scale from the Wall.
3. Check the measured sizes against canon (`world/canon/CONSTRAINTS.md`).

## Input data

- Map M1, 9050 × 6000 px.
- Script `scripts/land_mask.py`; parameters at the top of the script.
- Control points of the Wall, read on M1 by eye at full resolution: Westwatch (1568, 1015) and Eastwatch (1933, 962),
  where the Wall meets the Bay of Ice and the Bay of Seals. Reading error about ±10 px (±3 %).
- Run: `python3 land_mask.py --map <M1> --out <dir>` (about 20 s).

## Method and iterations

| Run | Change | Result |
|---|---|---|
| a | Sea = blue hue, saturation ≥ 0.18 | Forests and mountain shading (blue-grey) became sea; sea names stayed as land specks. |
| b | Saturation ≥ 0.30 | Forests and mountains fixed, but the pale shallow water along the coasts became land and joined Westeros to Essos. |
| c | Pale water (saturation ≥ 0.10, value ≥ 0.62) added only where it connects to confident sea; lettering removed by shape (strokes thinner than ~13 px) and colour (cream) | Coasts correct; ice floes partly land. |
| d | Ice floes = components lying wholly in the top 900 px; the Wall (a blue band of ice) forced to land | Floes removed. |
| e | Lettering test repeated after the morphological opening | Almost all sea names removed. Final run. |

Colour statistics that set the thresholds (median of sample windows): open sea saturation 0.35–0.54; pale coastal
water 0.36–0.40 with value up to 0.87; forest and mountain shading 0.19–0.27.

## Results (run e)

| Measure | Map M1 | Canon | Agreement |
|---|---|---|---|
| Wall | 368.8 px → **1.301 km/px** | 480 km (100 leagues, grade B) | defines the scale |
| Westeros, Wall to the south coast of Dorne | 4 000 km (2 490 mi) | ~3 000 mi for the whole continent, Land of Always Winter included (WE2) | consistent: the map is cut at the top and the north continues |
| Westeros, top edge of the map to Dorne | 5 150 km (3 200 mi) | ~3 000 mi (WE2) | +7 %, consistent |
| Westeros, widest east–west extent | 1 890 km (1 170 mi), at the latitude of the Iron Islands and the Fingers | ~900 mi (WE2) | **30 % wider** (see conclusions) |
| Westeros, area south of the Wall | 5.4 million km² | none stated | — |
| Westeros, area on the map | 7.0 million km² | "similar to South America" (17.8 million km², WE3) | not followed (deviation 2) |
| Essos, area on the map | 21.4 million km², 8 400 km east–west | "roughly Eurasia" (55 million km², E2) | consistent: Essos runs on beyond the east edge of the map |
| Whole map at the Wall's scale | 11 780 × 7 810 km | — | — |

The mask has 320 land components; the land is 38.6 % of the mapped area.

## Conclusions

- **Goal 1 reached.** The mask follows the canon coasts. Known defects:
  - islands smaller than ~150 px (~250 km²) are lost; this is below the 10 km grid (100 km² per cell) for most of them,
    but some named small isles may vanish;
  - the Thousand Islands and the Land of Always Winter keep some noise (blue shading read as water); real lakes of the
    Land of Always Winter are mixed with it;
  - inland lakes are treated like the sea when they are dark blue, and as land when small; lakes need their own layer.
- **Goal 2 reached.** 1.301 km per px at the Wall, ±3 %.
- **Goal 3 reached, with one open question.** North–south lengths and the size of Essos agree with canon. The width of
  Westeros measured on the map is ~1 170 mi, not ~900 mi. The 900 mi figure is A Wiki of Ice and Fire's own estimate
  from the Wall scale, not a statement in the books; the map is the canon source for shapes (`CONSTRAINTS.md`), so the
  map's width stands.
- **Main limitation: the map is flat.** The Wall gives the scale only where the Wall is. The map covers ~7 800 km north
  to south (about 70° of latitude on an Earth-sized planet), so on a globe its east–west scale must change with
  latitude, and distances far from the Wall (Essos, Sothoryos) are uncertain until the map's projection is decided.
  This is the next experiment: placing the map on the sphere.
