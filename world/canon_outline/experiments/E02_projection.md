# E02. Placing the flat map on the sphere

## Goal

Find how map M1 (`research/SOURCES.md`) can lie on an Earth-sized sphere (`world/planet/PLANET.md`) so that, as far
as possible:

1. the Wall keeps its canon length (480 km, the scale found in E01);
2. the outlines keep their canon shapes (owner's decision, `world/canon/CONSTRAINTS.md` section 7);
3. north stays "up" on the map, as the books' directions assume;
4. the canon climates fall in the latitude bands where the same climates occur on Earth (Rule 1).

A flat map 11 800 × 7 800 km cannot keep scale, shape and direction at once on a sphere (a sphere cannot be flattened
without distortion), so the experiment measures the trade-off.

## Input data

- The land mask from E01 (run e).
- Script `scripts/projection_study.py`. Run: `python3 projection_study.py --mask <land_mask.png> --out <dir>`
  (about 1 min). Placement scanned for the Wall's latitude 50–70° N in steps of 0.5°.
- Canon climate regions, boxes read on M1, and their Earth latitude bands:

| Region | Canon | Earth analogue and band | Source of the band |
|---|---|---|---|
| Haunted forest, north of the Wall | forest (WE7) | boreal forest, 50–70° N | [EBSCO, Boreal forests and taiga](https://ebsco.com/research-starters/earth-and-atmospheric-sciences/boreal-forests-and-taiga); Arctic treeline 60–70° N |
| Dorne | the only desert of Westeros, hottest region (WE9) | subtropical hot desert, 15–35° | [Britannica](https://www.britannica.com/science/tropical-and-subtropical-desert-climate) (most at 15–30°); [Physical Geology 2e](https://pressbooks.ccconline.org/accphysicalgeology/?p=1320) (15–35°) |
| Red waste | desert (E6–E7) | as Dorne, 15–35° | as above |
| Dothraki sea | grassland plain (E6) | temperate grassland, 30–60° | [Gervais, Midlatitude biomes](https://digfir-published.macmillanusa.com/gervais1e/gervais1e_ch08_4.html); [EBSCO, Temperate grasslands](https://www.ebsco.com/research-starters/earth-and-atmospheric-sciences/temperate-grasslands) |
| North coast of Sothoryos | jungle (S2) | tropical rainforest, 0–15° | [Wikipedia, Tropical rainforest climate](https://en.wikipedia.org/wiki/Tropical_rainforest_climate) |
| Summer Isles | tropical archipelago (I1) | tropics, 0–23.4° | Tropic of Cancer |
| Land of Always Winter | permanent ice (WE7) | ice sheet, ≥ 60° N | Greenland, 60–83° N |

## Ways of reading the map

The Wall fixes the scale in every case.

| | Projection | What it keeps | What it gives up |
|---|---|---|---|
| A | cylindrical equidistant, true at the Wall | north-up; north–south distances; parallels horizontal | east–west scale away from the Wall |
| M | Mercator, true at the Wall | north-up; shapes (conformal) | sizes change with latitude in every direction |
| B | equidistant conic, centre of the map | sizes and shapes within ~4 % | north tilts strongly towards the edges |
| C | equidistant conic, centred on the Wall | Westeros' shape and north-up | Essos tilts; Essos' rows bend to other latitudes |
| S | sinusoidal, centred on the Wall | north–south distances; **all** east–west distances; **all** areas; parallels horizontal | meridians lean (shear) away from the centre |
| S2 | sinusoidal, centre of the map | as S | as S, centred on Essos |

## Results

Best Wall latitude for each reading (the one where the most land of every climate region lies in its band):

| | Wall lat | Regions in band (worst) | Scale error, Westeros (95 %) | North tilt, Westeros (max) | Latitudes of the map | Longitudes spanned | Westeros area |
|---|---|---|---|---|---|---|---|
| A | 63.5° | **100 %** | 49 % (east–west) | 0° | 6–74° N | 231° | 9.7 M km² |
| M | 50° | 0 % | 35 % | 0° | 30° S–59° N | 161° | 11.2 M km² |
| B | 57.5° | 17 % | 3.5 % | 64° | 0–67° N | 166° | 7.2 M km² |
| C | 68.5° | 55 % | 17 % | 33° | 22° S–79° N | 149° | 7.5 M km² |
| **S** | **63.5°** | **100 %** | **7.6 %** | 41° (far north only) | 6–74° N | 179° | 7.0 M km² |
| S2 | 63.5° | 100 % | 127 % | 71° | 6–74° N | 248° | 7.0 M km² |

With A and S all seven regions fit their bands for a Wall latitude of **63.5–64° N** (≥ 90 % for 62–66.5° N). There,
Winterfell lies at ~56° N, Dorne at 28–32° N, the Dothraki sea at 30–44° N, the red waste at 23–32° N, the north
coast of Sothoryos at 6–14° N, the Land of Always Winter at 65–74° N.

How much S (centred on the Wall) bends each part of the map (tilt of local north; scale along the drawn vertical):

| Part | Tilt, median / 95 % | Scale, 95 % |
|---|---|---|
| Westeros south of the Wall | 3° / 11° | 1.7 % |
| Westeros beyond the Wall | 13° / 31° | 17 % |
| Western Essos (Free Cities) | 14° / 25° | 10 % |
| Central Essos (to the Bone Mountains) | 28° / 42° | 35 % |
| Eastern Essos (Jogos Nhai to Asshai) | 35° / 51° | 59 % |
| Sothoryos | 6° / 10° | 1.6 % |

## Conclusions

1. **The map was drawn with horizontal parallels.** Only the readings that keep the map's rows as parallels (A, S)
   put every canon climate in its Earth band; the conic readings (B, C) bend the rows of Essos to other latitudes
   (the red waste and the Dothraki sea fall out of their bands), and Mercator (M) spreads the map over 90° of
   latitude. Canon climate zonation is therefore the strongest evidence of how the map lies on the globe.
2. **The Wall lies at 63.5–64° N.** Westeros then runs from ~74° N (the map's top edge) to ~28° N.
3. **A and S differ in what they sacrifice.** A keeps north-up everywhere but makes the south of Westeros up to twice
   as wide as drawn (east–west scale at 30° N is cos 30° / cos 63.5° ≈ 1.95), Westeros' area grows from 7.0 to
   9.7 million km² and the known world spans 231° of longitude. S keeps every east–west distance and every area as
   drawn and keeps Westeros' shape (south of the Wall within 1.7 % and 11°), at the cost of shearing eastern Essos:
   there local north leans 30–50° from the map's vertical.
4. **Recommendation: S, centred on the Wall.** It is the only reading that keeps the Wall, the outlines of Westeros,
   all areas and the canon climates together; what it gives up (directions in eastern Essos) is the least-known part
   of canon. This needs the project owner's decision, because each reading gives up a different canon property.
5. **Not tested here:** whether the Haunted forest gets a warm enough summer at 63.5–66.5° N (Earth's treeline varies
   from ~58° N on cold coasts to ~70° N inland). The climate layer must check it; if it fails, the Wall moves south
   within 62–64° N, or the forest becomes a canon deviation.

## Decision

The project owner chose **S**, centred on the Wall. The Wall's midpoint is placed at **63.75° N**, the middle of the
range that fits every climate region. Recorded in `world/canon/CONSTRAINTS.md` (section 7) and, for the directions
it changes, in `CANON_DEVIATIONS.md` (entry 3).
