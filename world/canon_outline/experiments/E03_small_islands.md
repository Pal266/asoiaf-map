# E03. Small islands and map labels

## Goal

1. Restore the islands that E01 lost: E01 erased every land body thinner than ~9 px or smaller than 150 px
   (~250 km²) with its opening and size filter, together with the lettering it was meant to remove.
2. Keep the islands on the global level, so that regional levels can build them (Rule 6: a lower level is built only
   from the level above, so an island missing here would be missing everywhere below).
3. Store the islands too small for the grid to draw well as polygons.

## Input data

- Map M1, the classification of E01 before its opening (`raw_land`) and its final mask.
- Script `scripts/land_mask.py`, steps 6–7 (run as in `GENERATION.md`).

## Method and iterations

Large land masses (≥ 20 000 px, ~34 000 km²) keep the E01 outline. Every smaller body is rebuilt.

| Run | Rule | Result |
|---|---|---|
| a | Candidates = E01 raw bodies that its opening erased; inspected on contact sheets | 5 951 candidates; most were sea speckle (dark or weakly saturated sea), lettering and fold lines; real islands were few and small (18–150 px) |
| b | Island = a bright core whose surrounding ring is darker than the sea (the coastline stroke) | Dark map labels (island and place names: dark ink with a light halo) passed as islands |
| c | Core must have land-fill colours | Labels still passed (ink anti-aliasing looks tan) |
| d | Ring must be olive-brown (outlines are olive, label ink is purple-maroon) | Labels gone, but most islands lost too (outline colour varies) |
| e | Remove label ink (dark purple-maroon) and halos (light, weakly blue) first, then keep a raw body if its interior holds land-fill colours | Islands of every group found; cream sea names remained |
| f | Remove weakly saturated sea (blue, saturation ≥ 0.2) from candidates; letters = thin bodies that are mostly pale | Speckle gone; cream letters gone, but some pale grey islets lost |
| g | A thin pale body with an edge darker than the sea is kept; pink letters added to "pale"; purple counted as land only in the Smoking Sea (elsewhere mauve is lettering) | Lettering gone, grey islets back |
| h | Weak-sea saturation 0.25 (grey-blue rock shading of islands reaches ~0.2); Old Wyk added from a box (it is pale with blue-grey shading and no darker edge on the navy sea) | Final run |

Colour evidence (pixel samples on M1): label ink hue 266–347°, saturation 0.1–0.3, value 0.25–0.33; island outline
hue 44–53°, value 0.31–0.40; island fill khaki/green (hue 40–80°, saturation 0.3–0.45) or grey rock (saturation
< 0.25); Valyrian land dark purple; cream sea names hue 35–50°, saturation 0.15–0.3, value 0.75–0.9, and greenish
(hue 100–180°, saturation ~0.13) where blended into light sea.

## Results (run h)

- **286 islands smaller than 2 500 km²** traced as polygons, 161 000 km² in all; 23 are smaller than one grid cell
  (100 km²); the smallest is 44 km².
- Against E01: 87 600 km² of land added (islands) and 163 500 km² removed (map labels and sea speckle that E01 had
  kept as land: the names of the Iron Islands, Bear Island, Skane, Lorath Bay, the Basilisk Isles, the copyright line
  and others). Land components: 359 (E01: 320).
- Checked by eye on the Iron Islands (Great Wyk, Old Wyk, Pyke, Saltcliffe, Harlaw, Orkmont, Blacktyde), the
  Stepstones, the Summer Isles, the Thousand Islands, the Basilisk Isles, Valyria and the Smoking Sea, Skagos and
  Skane, Ibben, the Isle of Cedars, Lorath, Bear Island.
- On the grid (8 × 8 samples per cell, about one map pixel apart): land area matches the map to < 0.001 %.

## Conclusions

- **Goal 1 reached** for islands with a land-coloured fill. Still lost: islets drawn as dark specks without fill (many
  in the Thousand Islands and along some coasts, mostly under ~40 km²) and any label that lies over a large land
  mass's coast (it stays part of that coast, as in E01).
- **Goal 2 reached:** the islands are part of the global land mask and of the grid layer.
- **Goal 3 reached:** `islands.geojson` holds each island's polygon (longitude/latitude), area, centroid and centroid
  cell, so regional levels can rebuild islands smaller than the grid resolves.
- The classifier needed one manual box (Old Wyk). Further named islands found missing later are added to
  `MANUAL_ISLANDS` in the same way.
