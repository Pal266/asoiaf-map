# E12 — Hydrology: runoff, rivers and lakes

Status: **done** — run k is the layer (decisions 9 and 10, canon deviation 7).

## Goal

Annual runoff, river discharge and lakes for every cell, from the climate with the ice (E09–E11) and the ice
surface. Canon checks:
- the Rhoyne drains much of western Essos and reaches the Summer Sea at Volantis (E8);
- several large lakes separate the Land of Always Winter from the haunted forest (WE7);
- the Shrinking Sea and the Bleeding Sea are lakes;
- the rivers drawn on the official map.

## Model (`scripts/flow.py`)

- **Potential evapotranspiration:** Thornthwaite (1948) from monthly temperature and day length, × `pet_factor`.
- **Runoff:** the Budyko curve in Fu's form, E/P = 1 + φ − (1 + φ^w)^(1/w), φ = PET/P, w = 2.6 (Zhang et al.
  2004). On ice sheets all precipitation runs off.
- **Routing:**
  - priority-flood filling from the sea (Barnes et al. 2014);
  - steepest descent over 8 neighbours on the filled surface;
  - accumulation downstream, never below zero.
- **Lakes:** depressions of the unfilled surface.
  - If the inflow exceeds the open-water net evaporation of the whole depression, it overflows. Only water deeper
    than 50 m stays lake: the shallow margins of an overflowing basin fill with sediment.
  - Otherwise it is a terminal lake, shrunk to its lowest cells until their net evaporation equals the inflow; its
    river ends there.

## Earth test (`scripts/earth_hydro_validation.py`)

- **Inputs:** Earth climate from the E09 and E10 models at 10′ (as in the ice test) and ETOPO1.
- **Rivers:** 16 large rivers with gauge discharge ("List of rivers by discharge", Wikipedia).
- **Stream burning:** ETOPO1 at 10′ comes from a smoothed spherical-harmonic model, which loses the gorges that carry
  rivers. Without help the Congo drained into the Nile. So Natural Earth 1:50m river centrelines are burnt into it
  (1000 m), and lakes are looked for on the unburnt surface.

| Run | Change | log RMSE of Q | drainage area log RMSE | runoff to sea km³/yr | Notes |
|---|---|---|---|---|---|
| a | first | 1.56 | 1.06 | 48 800 | the Congo is 177 000 km² and drains into the Nile |
| b | stream burning | 1.18 | 0.65 | 50 700 | 1 M km² of terminal "lakes" along burnt channels |
| c | lakes on the unburnt surface | 0.71 | 0.36 | 50 750 | full lakes of ~1 M km² in smoothed flat basins (Congo, west Siberia) evaporate the rivers |
| d | overflowing lakes only deeper than 50 m | 0.67 | 0.36 | 52 200 | |
| e | PET × 1.5 | 0.83 | 0.36 | **40 300** | global runoff to the sea as on Earth (~37 300 km³/yr of rivers, Dai & Trenberth 2002, plus ice-sheet melt) |

**The river errors come from the precipitation, not the water balance.** Basin means of model precipitation
against the observed runoff (mm/yr):
- Mississippi: P 217, runoff 210 observed (the eastern United States are ~3× too dry in E10);
- Irrawaddy: P 689, runoff 1185 observed (the monsoon is too weak);
- Nile: P 1260, runoff 30 observed (too wet; and the Sudd's evaporation is not modelled);
- Lena, Mackenzie and Yangtze, where precipitation is right, are within 10–60 %.

PET × 1.5 is the one fitted parameter: Thornthwaite's formula underestimates potential evaporation in arid and
tropical climates. It is set by Earth's total runoff to the sea.

## World, first run (a)

- Land runoff 298 mm/yr on average; 38 300 km³/yr reaches the sea.
- Lakes: 5.0 M km² overflowing plus 8.4 M km² terminal. Earth has ~5 M km² of lakes of all sizes, ~1.5 of them
  large. A terminal inland sea of 5.3 M km² fills a closed basin on the south-western continent; its inflow is
  200 000 m³/s.
- Canon checks fail:
  - Rhoyne at Volantis 17 m³/s (basin 81 000 km²);
  - Shrinking Sea and Bleeding Sea 1–2 % lake (P ~100–140 mm/yr, no inflow);
  - Land of Always Winter lakes 1 500 km².
- The relief (E08) has no river incision. Its drainage follows tectonic relief and noise, and does not reproduce the
  rivers on the official map: in the Riverlands the main flow goes west, while the canon Trident flows east.

**Owner decision 1 (2026-11-22): canon rivers into the relief** (CONSTRAINTS decision 9).

## Canon rivers (`scripts/canon_rivers.py`, `data/canon_rivers.json`)

- **26 rivers traced by eye** on gridded views of M1, as polylines in map pixels:
  - Westeros: the Trident and its three forks, the Blackwater Rush with the God's Eye outflow, the White Knife,
    the Weeping Water, the Mander with four tributaries, the Honeywine, the Wendwater, the Greenblood with the
    Scourge;
  - Essos: the Rhoyne with nine tributaries.
- **Snapping:** each vertex moves to the darkest pixel within 12 px, then each segment follows the least-cost path
  through the map's dark pixels inside an 18 px corridor (`check.png`).
- **Corrections:** in the first pass the Mander and the Wendwater followed the dotted Roseroad and the Honeywine a
  neighbouring line. They were re-traced on closer views. The trace first named "Last River" is the Weeping Water.
- **Carving** (`hydrology.carve_canon`): canon_burn_m (500 m) below the surface, falling ≥ 0.5 m per cell from
  source to mouth. A floor of 0.5 m × (cells to the sea along the network) keeps each river above sea level and above
  the river it joins. A mouth that stops short of the model's coast is continued to the nearest sea cell. Lakes are
  still found on the uncarved surface.

| Run | Change | Rhoyne at Volantis | Notes |
|---|---|---|---|
| b | burnt 500 m | 18 m³/s | the network is followed, but the burnt channels drain out at whichever end is lowest |
| c | falling from source to mouth | 11 | the Rhoyne drains up the Darkwater: that tributary's end lay below the Rhoyne's junction cell |
| d | backward pass above junctions | 21 | channels raised above the burn level leak sideways |
| e | floor by distance to the sea | 11 | the Rhoyne's traced mouth ends on land at 341 m; the pit spills north through the Qhoyne's low source |
| **f** | **mouths continued to the coast** | **1169 m³/s, basin 1.36 M km²** | every canon river drains along its traced course (84–100 % of its cells) |
| g | marsh and salt flat classes | — | below |

**Marshes and salt flats** (run g):
- the sediment-filled part of an overflowing depression (1–50 m deep) is marsh or floodplain: 3.6 M km², 23 % of
  the Neck (canon: wetland);
- the floor of a terminal depression up to 20 m above its lake is salt flat: 2.4 M km². Counting every cell below
  the spill point gave 15 M km², most of it alluvial plain.

## World (run g)

| River | Discharge m³/s | River | Discharge m³/s |
|---|---|---|---|
| Trident (with the forks) | 3 150 | Rhoyne at Volantis | 1 170 |
| White Knife | 3 990 | Mander | 230 |
| Weeping Water | 1 080 | Blackwater Rush | 100 |
| Greenblood | 8 | Honeywine | 1 |

For scale: the Danube 6 500, the Rhine 3 200, the Guadalquivir ~160.

- Land runoff 298 mm/yr; 38 300 km³/yr to the sea.
- Lakes 5.0 M km² overflowing plus 8.4 terminal, mostly outside the known world: a 5.3 M km² inland sea in a closed
  basin of the south-western continent, and lakes of 1.6 and 1.2 M km² in deep basins of Sothoryos.
- The North has two lakes of ~20 000 km² near 60° N.
- **The south of Westeros is dry.** The Reach gets ~300–400 mm/yr against high potential evaporation (E10: a
  Mediterranean regime, King's Landing 268 mm), so the Mander carries 230 m³/s from 600 000 km² (12 mm/yr of runoff).
- **Canon lakes:**
  - the Shrinking Sea is 1 % lake and 12 % salt flat;
  - the Bleeding Sea is 2 % lake;
  - their basins get 94–141 mm/yr and almost no inflow.
- The Land of Always Winter has 1 500 km² of lakes; canon (WE7) has several large ones.

**Owner decisions 2 and 3 (2026-11-22): the canon lakes take their water from the map; revisit the precipitation.**

## Canon lakes (run h)

Envelopes around four lakes drawn on M1 (`data/canon_rivers.json`, `lakes`). The water inside is blue, pale grey-blue
(Dagger Lake) or red (the Bleeding Sea), and it becomes lake whatever the depressions say (`forced_lake`). The
open-water balance then decides whether a lake overflows.

| Lake | Area km² | Outflow m³/s | |
|---|---|---|---|
| God's Eye | 13 400 | 39 | overflows to the Blackwater Rush |
| Dagger Lake | 2 300 | 1 130 | on the Rhoyne |
| Shrinking Sea | 7 900 (two pools) | 0 | terminal, 11 % of its envelope salt flat: the map's dried brown bed around two pools |
| Bleeding Sea | 96 400 | 0 | terminal: its evaporation exceeds its inflow |

## Revisiting the precipitation (run p10, reverted)

The E10 model was recalibrated with the Earth river test added to its objective: discharge of 16 rivers, and
runoff to the sea of 39 000 km³/yr. 14 parameters, including `pet_factor`, over 320 evaluations:
- Earth improved a little: stations log RMSE 0.388 → 0.376 (61 % within a factor 2); rivers 0.71 → 0.63
  (Mississippi −1.17 → −0.94); arid stations worse (+0.25 → +0.34);
- the Reach did not change: King's Landing 268 → 269 mm, Oldtown 325 → 331 mm. Its dryness is robust: a west coast at
  30–40° N on a Pacific-sized ocean, like California, Chile or Morocco;
- the North got wetter (Winterfell 689 → 789 mm), and the ice layer rerun grew a northern ice sheet again
  (14.6 → 16.3 M km² and growing): the bistable northern ice (E11) flipped.

So the recalibration was reverted; parameters stay as committed.

## Lakes in noisy relief (runs i–k)

Small closed pits of the relief became thousands of small terminal lakes and salt flats in dry interiors. Depressions
under 1000 km² or less than 10 m deep at their deepest point are now taken as filled by alluvium (run k). Using a
10 m threshold on the depression's extent instead (run j) split depressions into parts and was rejected.

## World (run k)

| River | m³/s | River | m³/s |
|---|---|---|---|
| White Knife | 3 990 | Rhoyne at Volantis | 1 160 (basin 1.36 M km²) |
| Trident | 3 120 | Weeping Water | 1 080 |
| Mander | 240 | Blackwater Rush | 94 |
| Greenblood | 8 | Honeywine | 1 |

- Land runoff 298 mm/yr; 38 300 km³/yr to the sea.
- Lakes: overflowing 5.0 M km², terminal 8.4, marsh 3.4, salt flat 2.3. Most of the lake area lies outside the known
  world: the south-western inland sea (5.3 M km²) and two deep basins in Sothoryos.
- The Neck is 23 % marsh.

## Canon conflict (Rule 3), resolved

The Reach and the Crownlands have a dry Mediterranean climate in every version of the model. So the Mander (240 m³/s
from 600 000 km²) and the Blackwater Rush (~100 m³/s) are much smaller rivers than the map and the books suggest.

**Owner decision 4 (2026-11-22): canon deviation 7.** The dry Mediterranean south stands, with smaller rivers.

## Still open

- Seasonal discharge: the layer has annual values only. A monthly snow and soil store is the next step if the game
  needs floods and low water.
- Flow over large flats follows parallel lines (priority-flood artefact); marshes and salt flats mark most of them.
- The relief has no valleys: the canon rivers are carved only in the routing surface, not in `elevation_m`.
- Big lakes and an inland sea outside the known world come from closed basins of the generated relief.

## Sources

- Thornthwaite, C. W. (1948), Geogr. Rev. 38, 55 — potential evapotranspiration; Willmott et al. (1985) above 26.5 °C.
- Fu, B. (1981); Zhang, L. et al. (2004), Water Resour. Res. 40 — Budyko curve, w ≈ 2.6.
- Barnes, R., Lehman, C. & Mulla, D. (2014), Comput. Geosci. 62, 117 — priority flood.
- Saunders, W. (1999) — stream burning in DEM conditioning.
- Dai, A. & Trenberth, K. E. (2002), J. Hydrometeorol. 3, 660 — continental discharge, 37 300 km³/yr.
- "List of rivers by discharge", Wikipedia (read 2026-11-22).
- Natural Earth 1:50m rivers and lake centrelines (public domain).
- Map M1, *The Lands of Ice and Fire* (2012) — canon rivers and lakes (not stored).

