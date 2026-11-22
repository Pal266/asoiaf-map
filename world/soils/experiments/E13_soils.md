# E13 — Soils: a climate and terrain key to the WRB Reference Soil Groups

Status: **done** — run c is the layer.

## Goal

The dominant soil of every cell, as a WRB Reference Soil Group, and a natural fertility index for the game. It is built
from the layers above: climate with the ice, relief, ice, hydrology and canon volcanic relief.

## Model (`scripts/soilkey.py`)

Soils are named by the factor that most clearly conditions them. The FAO lecture notes (World Soil Resources Report
94, 2001) group the Reference Soil Groups in sets that way. The key follows the sets in order:

1. **Water, ice, organic and salty:**
   - lakes → water; ice sheets → ice;
   - marshes in cool, moist climates and wet flats of the boreal lowlands → Histosols;
   - salt flats → Solonchaks.
2. **Permafrost:** mean annual temperature ≤ −2 °C → Cryosols.
3. **Topography and parent material:**
   - strong local relief, or moderate relief in arid land → Leptosols;
   - canon volcanic relief → Andosols;
   - rivers ≥ 500 m³/s → Fluvisols;
   - other marshes → Gleysols.
4. **Deserts** (Thornthwaite arid and < 250 mm): Gypsisols in the driest closed basins, Arenosols on warm low
   plains, Leptosols below 50 mm, Calcisols otherwise.
5. **Humid tropics, old stable surfaces:** coldest month ≥ 18 °C, moist, low relief, dry season ≤ 5 months →
   Ferralsols.
6. **Seasonal clay plains:** basins and lowlands with ≥ 4 dry months → Vertisols.
7. **Sandy plains of the dry tropics with summer rains** (Sahel, Kalahari) → Arenosols. With winter rains a dry
   plain is Mediterranean instead.
8. **Steppes with cold winters**, by moisture class:
   - semi-arid → Kastanozems (FAO: 200–400 mm);
   - dry subhumid → Chernozems (FAO: precipitation ≈ evaporation);
   - moist subhumid with warm summers → Phaeozems.
9. **Wet tropics:** Umbrisols above 2000 m, Acrisols if moist, Lixisols if drier.
10. **Subtropics and Mediterranean:** Calcisols (semi-arid), Acrisols (humid), Cambisols (hilly), Luvisols.
11. **Cool temperate and boreal:** Podzols, Albeluvisols (cold continental), Umbrisols (perhumid oceanic),
    Cambisols (hilly), Kastanozems (dry), Luvisols, Cambisols.

**Moisture classes** are Thornthwaite's (1948), from Im = 100 (P / PET − 1), with plain Thornthwaite PET:
- arid < −40; semi-arid −40 to −20; dry subhumid −20 to 0;
- moist subhumid 0 to 20; humid 20 to 100; perhumid > 100.

The ×1.5 factor of the hydrology layer is a runoff calibration and is not used here.

**Local relief** is the height range over ~50 km. The thresholds are set on the Earth test's grid (3 × 3 cells of
10′ ETOPO1). The world's 10 km relief (5 × 5 cells) has a different spread: median 115 against 100 m, 10th
percentile 41 against 20 m, 97th 718 against 960 m. So it is mapped onto Earth's distribution by rank (quantile
mapping, `earth_relief_quantiles_m`).

**Fertility** (`soils.json`, `fertility`) is a 0–1 ranking for the game, from the agricultural use described for each
group in the lecture notes:
- Chernozems 1.0, Phaeozems 0.9, Fluvisols 0.85, Luvisols and Andosols 0.75;
- Podzols 0.25, Arenosols 0.15, Leptosols 0.1, Solonchaks 0.05.

This is a judgement, not a measurement.

## Earth test (`scripts/earth_soil_validation.py`)

- **Inputs:** Earth climate from the E09 and E10 models at 10′, ETOPO1, today's ice, and the Earth hydrology test
  (marshes, salt flats, discharge). No volcanic map, so Andosols are not tested.
- **46 sites** (`data/earth_soil_sites.json`) in regions that the lecture notes name as major areas of a group: the
  Brazilian and African shields for Ferralsols, the Ukrainian steppe for Chernozems, Mongolia for Kastanozems,
  Sudan and the Deccan for Vertisols, Yakutia for Cryosols and so on. A site counts if an accepted group occurs at
  its cell or a neighbour.
- **Areas** of each group against the FAO estimates in the same notes.

| Run | Change | Sites | Area log RMSE | Notes |
|---|---|---|---|---|
| a | first key, aridity index with PET × 1.5 | 41 % | 0.54 | Vertisols 1914 against 335 Mha; Cambisols 175 against 1500 |
| b | basins and lowlands for Vertisols, relief classes, broader Ferralsols | 48 % | 0.40 | |
| c0 | Ferralsols before Vertisols | 50 % | 0.40 | |
| d | Thornthwaite moisture classes with plain PET; FAO precipitation ranges for steppes and deserts | 52 % | 0.28 | |
| **e** | relief thresholds adjusted; peat flats only below 50 m of relief | **52 %** | **0.26** | chosen |

Areas, model against FAO (Mha):

| Group | Model | FAO | Group | Model | FAO |
|---|---|---|---|---|---|
| Cryosols | 1863 | 1800 | Ferralsols | 767 | 750 |
| Leptosols | 1114 | 1655 | Acrisols | 680 | 1000 |
| Cambisols | 467 | 1500 | Lixisols | 411 | 435 |
| Calcisols | 1460 | 1000 | Luvisols | 408 | 550 |
| Arenosols | 1991 | 900 | Podzols | 354 | 485 |
| Gleysols | 272 | 720 | Histosols | 393 | 350 |
| Vertisols | 678 | 335 | Kastanozems | 817 | 465 |
| Chernozems | 456 | 230 | Phaeozems | 73 | 190 |

**Most misses come from the Earth climate model, not from the key.** At the sites the model gives:
- southeastern USA 354 mm, Great Plains 103 mm, Iowa 252 mm, Texas 121 mm (observed ~1200, 450, 850, 900);
- the Sahel and the Sudan plains ~1000 mm (observed ~400);
- northwest Russia an aridity that puts it in the steppe.

These are the E10 faults (E10, E12). Where the model climate is right, the key gives the textbook soil:
Ukrainian and Russian steppe Chernozems, Mongolia and southern Russia Kastanozems, Scandinavia and Canada Podzols,
Yakutia, Alaska and northern Canada Cryosols, Congo and Brazil Ferralsols, Kalahari and Australia Arenosols.

## World (run c)

| Region | Main groups (share of land) | Mean fertility |
|---|---|---|
| Westeros south of the Wall | Calcisols 19 %, Podzols 14 %, Luvisols 12 %, Histosols 10 %, Cambisols 10 %, Kastanozems 7 %, Cryosols 6 %, Arenosols 6 % | 0.38 |
| Beyond the Wall | Cryosols 94 %, Histosols 5 % | 0.11 |
| Western Essos | Calcisols 29 %, Luvisols 16 %, Kastanozems 16 %, Leptosols 13 %, Cambisols 10 % | 0.47 |

- **The North:** Podzols with peat on the low wet ground, Cryosols in the colder uplands.
- **The Riverlands:** Luvisols and Cambisols.
- **The Neck:** Gleysols (marsh).
- **The Reach and the Crownlands:** Kastanozems and Calcisols: the dry Mediterranean south (canon deviation 7).
- **Dorne:** Arenosols and Calcisols, with Leptosols in the Red Mountains.
- **Valyria:** Andosols (volcanic); the canon volcanic features of the Fourteen Flames give 0.33 M km² of Andosols.
- **Elsewhere:** the Rhoyne and the other large rivers have Fluvisol belts; beyond the Wall the permafrost of the
  colder climate (E11) makes Cryosols.

## Still open

- ~~Parent material is known only where it is volcanic.~~ Done: the lithology layer (E14), see the addendum.
- No Regosols, Planosols, Solonetz, Plinthosols, Alisols, Durisols or Anthrosols.
- A single dominant group per 10 km cell; no soil associations or properties (texture, depth, pH).

## Addendum: parent material from the lithology (run d)

The key takes `rock`, `deposit` and `bedrock` from the lithology layer (E14) when it is given:

- **Andosols:** young volcanics and pyroclastics in moist climates (FAO: humid volcanic regions), besides the canon
  volcanic relief.
- **Arenosols:** aeolian sand in deserts and semi-arid lands; the relief-based sand rules apply only without
  lithology.
- **Nitisols** (new group, fertility 0.8): moist tropics on basalts and gabbros, on moderate relief, before
  Ferralsols.
- **Vertisols:** warm, seasonally dry, low clay plains: basalts, alluvium, closed basins and mixed sediments
  (shales, marls).
- **Leptosols:** moist limestone uplands (rendzinas).
- **Podzols and Luvisols:** moist sand and moist loess respectively.

| Group (M ha) | run c | run d |
|---|---|---|
| Arenosols | 1470 | 895 |
| Vertisols | 294 | 74 |
| Andosols | 33 | 176 |
| Nitisols | 0 | 21 |
| Leptosols | 1380 | 1518 |
| Luvisols | 609 | 688 |
| Calcisols | 3231 | 3501 |

Arenosols fall because only true sand seas count now; dry plains without sand become Calcisols. Vertisols shrink
to the clay plains. The Earth test is unchanged (52 %, log RMSE 0.26): it runs without lithology.

| Region (run d) | Main groups | Mean fertility |
|---|---|---|
| Westeros south of the Wall | Calcisols 24 %, Podzols 13 %, Luvisols 13 %, Histosols 10 %, Cambisols 9 %, Kastanozems 7 %, Leptosols 7 % | 0.39 |
| Beyond the Wall | Cryosols 94 %, Histosols 5 % | 0.11 |
| Western Essos | Calcisols 28 %, Luvisols 17 %, Kastanozems 16 %, Leptosols 13 %, Cambisols 10 %, Arenosols 6 % | 0.47 |

## Sources

- FAO (2001), Lecture notes on the major soils of the world, World Soil Resources Report 94 (Driessen, Deckers,
  Spaargaren, Nachtergaele, eds.), fao.org/4/y1899e — environments, regional distribution and areas of the groups.
- IUSS Working Group WRB, World Reference Base for Soil Resources, 4th edition (2022) — the groups.
- Thornthwaite, C. W. (1948), Geogr. Rev. 38, 55 — potential evapotranspiration and moisture index.
