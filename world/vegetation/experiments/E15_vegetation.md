# E15 — Vegetation

**Goal:** the potential natural vegetation (the vegetation without people) of every 10 km cell as a biome, with
its net primary productivity. Built from the climate, soils, hydrology, relief and ice layers. It is the starting
landscape of the game (forest to clear, grassland to graze) and, later, the base of the biological resources.

## Research

- **The BIOME family of models** (Prentice et al. 1992; Haxeltine & Prentice 1996; Kaplan et al. 2003). Plant
  functional types are admitted by climatic limits:
  - winter cold: the mean of the coldest month, Tc; tropical trees need Tc ≥ 15.5 °C (BIOME1);
  - summer heat: growing degree-days above 5 °C (GDD5) and 0 °C (GDD0); boreal trees need ~350 GDD5;
  - moisture: α = actual / potential evapotranspiration from a 150 mm soil bucket.

  The biome is what the admitted types make together. The mega-biome limits of BIOME4, as tabulated by Dallmeyer
  et al. (2019), are:
  - tropical forest: Tc > 15.5 °C;
  - warm-temperate forest: GDD5 > 3000;
  - temperate forest: GDD5 900–3000;
  - boreal forest: GDD5 ≤ 900;
  - savanna: GDD5 > 1200 and Tc > 10 °C;
  - tundra: GDD0 < 800.
- **Mediterranean climate:** Köppen's Cs (the driest summer month < 40 mm and < 1/3 of the wettest winter month).
- **Savanna:** keeps an open canopy where the dry season is long, with fire and drought, even at moderate rainfall.
- **Productivity:** the Miami model (Lieth 1975), NPP = min(f(T), f(P)), in grams of dry matter. Global terrestrial
  NPP is ~56 Pg C/yr (Field et al. 1998), ~120 Pg of dry matter.
- **Earth reference:** the 14 biomes of Olson et al. (2001, WWF), with areas tabulated by Bertzky et al. (2013, IUCN):
  - deserts and xeric shrublands 28.0 M km²;
  - tropical grasslands and savannas 20.3;
  - tropical moist forests 19.9;
  - boreal forests 15.1;
  - temperate broadleaf and mixed forests 12.8;
  - tundra 11.6;
  - temperate grasslands 10.1;
  - montane grasslands 5.2;
  - temperate conifer forests 4.1;
  - Mediterranean 3.2;
  - tropical dry forests 3.0;
  - flooded grasslands 1.1;
  - tropical conifer forests 0.7;
  - mangroves 0.35.

## Method (scripts/vegkey.py)

1. **Bioclimate:**
   - Tc, Tw, MAT, GDD5 and GDD0 from the monthly temperature;
   - Thornthwaite potential evapotranspiration, as in the soil layer;
   - a Thornthwaite–Mather bucket of 150 mm, spun up three years, gives the actual evapotranspiration and α;
   - Köppen's Cs test and the winter half-year's share of the rain;
   - the number of dry months (P < PET/2).
2. **The key, in order:**
   - water and ice;
   - polar desert (GDD0 < 100);
   - tundra, or alpine below 55° latitude above 1500 m (GDD5 < 350);
   - desert (α < 0.15), xeric shrubland (α < 0.30);
   - tropics (Tc ≥ 12 °C):
     - savanna (α < 0.45 or ≥ 6 dry months);
     - tropical dry forest (α < 0.75);
     - tropical rainforest;
   - Mediterranean woodland: Cs, ≥ 60 % of the rain in the cool half-year, Tc ≥ 2 °C, α < 0.65;
   - boreal zone (GDD5 < 900, or Tc < −15 °C with GDD5 < 1500) where α ≥ 0.35:
     - forest-tundra (GDD5 < 600);
     - boreal forest;
   - temperate grassland: α < 0.5, or the dry boreal zone;
   - subtropical evergreen forest (Tc ≥ 5 °C, GDD5 ≥ 3000);
   - temperate conifer forest: above 800 m, or a cool, very wet coast (P ≥ 1500 mm, Tw < 19 °C);
   - temperate mixed forest (GDD5 < 1200 or Tc < −15 °C);
   - temperate broadleaf forest.
3. **High ground:** open land (tundra, grassland, xeric shrubland, savanna) above 2500 m below 55° latitude is
   alpine (WWF's montane grasslands).
4. **Water and soils:**
   - marshes of the hydrology and Histosols (peat) of the soil layer become wetlands;
   - where the climate carries forest with GDD5 ≥ 900, they become swamp forest instead;
   - rivers of ≥ 100 m³/s through desert, xeric shrubland, grassland or savanna are lined with riparian woodland;
   - salt flats stay bare.
5. **Productivity:** NPP from the Miami model, 0 on water, ice and salt flats.

## Earth test (scripts/earth_veg_validation.py)

The test uses the Earth inputs of the soil test: the model's own climate on Earth at 10′, ETOPO1, today's ice, the
Earth hydrology and the soil key's Histosols. So the thresholds absorb the model's errors on Earth in the way the
world run sees them.

It has two scores:
- **sites:** 78 textbook locations with the biome or biomes accepted there. 62 cover all biomes, and 16 lie along
  the arctic tree line. A site counts as right if an accepted biome is at its cell or one of the 8 neighbours.
- **areas:** WWF biome areas. Wetlands count as boreal or tundra below GDD5 900, swamp forests count as their
  climatic forest, and the other wetlands count as flooded grasslands. The score is the log10 RMSE over the biomes
  the key makes (no mangroves, no tropical conifer forests).

| Run | Change | Sites | Area log RMSE |
|---|---|---|---|
| a | BIOME1 starting limits (Tc 15.5, α 0.18 / 0.33 / 0.65, tropical moist 0.85) | 0.50 | 0.395 |
| b | Tc of the tropics 12 °C, moisture limits lowered by ~0.1, swamp forest, Mediterranean needs winter rain | 0.66 | 0.255 |
| c | savanna by dry-season length; cold-winter boreal; montane grassland | 0.65 | 0.230 |
| d | grid search of the tropical, forest and dryland limits (tree line left free: 650 GDD5) | 0.66 | 0.095 |
| e | temperate rainforest conifers; montane conifers from 800 m | 0.68 | 0.077 |
| f | + 16 tree-line sites (78 in all), tree line 650 | 0.667 | 0.078 |
| **g** | **tree line 350 GDD5 (BIOME1), forest-tundra to 600** | **0.641** | **0.142** |

- **Why the limits moved.** The model's climate on Earth has winters 5–9 °C too cold at the tropics' margins
  (Nagpur Tc 12 °C against ~21 °C, Brasília 14 against ~19, Chiang Mai 14 against ~21). And α comes from
  Thornthwaite's evapotranspiration, not BIOME1's Priestley–Taylor, so its limits were fitted again.
- **The tree line is not calibrated on Earth.** The model's arctic summers err in opposite directions, so model
  GDD5 does not separate forest from tundra:
  - too cold in Alaska and the Yukon: Fairbanks 413, Inuvik 323, Whitehorse 251 GDD5, all forest;
  - too warm around Hudson Bay and Baffin Island: Iqaluit 673, Baker Lake 556, Churchill 636 GDD5, all tundra.

  The area score prefers 650, the tree-line sites cannot decide, and the world's North would be tundra at 650
  even at Winterfell. So the forest limit is the physiological one of BIOME1, 350 GDD5.
- **Final run g:**
  - areas within a factor 1.3 of WWF for the tropical forests, temperate broadleaf and conifer forests,
    savannas, temperate grasslands, Mediterranean and deserts;
  - boreal forests 20.1 M km² against 15.1, flooded grasslands 1.5 against 1.1, tundra 5.3 against 11.6, montane
    grasslands 2.8 against 5.2;
  - global NPP 115 Pg of dry matter per year (~120).
- **Remaining misses** (28 of 78) follow the model's climate errors on Earth:
  - the wet Sahel (Kano, Niamey) is rainforest;
  - the too-dry US plains (Wichita, Cheyenne) are xeric shrubland;
  - Atlanta, Tokyo and Buenos Aires get too little rain;
  - three coastal sites fall on sea cells;
  - the arctic sites above.

The full test from the raw Earth data takes ~70 s and is deterministic.

## The climate of the North (CONSTRAINTS decision 11)

The first world run made the North of Westeros forest-tundra and bog and the Haunted Forest tundra. This goes
further than canon deviation 5 ("sparse taiga"). The owner had the climate of the North re-examined (E09,
addendum 2):
- the model's summers have no cold error at Earth's maritime subarctic stations;
- the cool summers come from the cold seas on both sides of the North and the polar sea's summer ice.

The owner kept the climate and extended deviation 5.

## World (run g)

| Region | Main biomes (share of land) | Mean NPP, g/m²/yr |
|---|---|---|
| Westeros south of the Wall | Mediterranean woodland 16 %, temperate broadleaf forest 16 %, xeric shrubland 11 %, wetland 10 %, forest-tundra 9 %, tundra 9 %, boreal forest 5 %, savanna 5 %, conifer forest 5 %, desert 5 % | 697 |
| Beyond the Wall | tundra 90 %, wetland 6 %, polar desert 4 % | 253 |
| Western Essos | temperate grassland 33 %, xeric shrubland 29 %, Mediterranean woodland 14 %, temperate broadleaf forest 12 % | 629 |

- **The North:**
  - 60–63.75° N: tundra 73 %, bog 22 %;
  - 56–60° N: bog 35 %, tundra 31 %, forest-tundra 30 %;
  - 52–56° N: bog 39 %, forest-tundra 38 %, taiga 15 %;
  - Winterfell stands among bog and forest-tundra (canon deviation 5).
- **The Neck:** forests with swamp forest on 23 % (the hydrology's marshes); closed forest begins here.
- **The Riverlands, the Vale and the Westerlands:** temperate broadleaf forest; the Mountains of the Moon carry
  conifer forest, taiga and alpine meadow.
- **The Reach and the Crownlands:** Mediterranean woodland (King's Landing), the dry south of canon deviation 7.
- **Dorne:** xeric shrubland and dry savanna, desert in the west; the Dornish desert is 54 % xeric shrubland, 34 %
  dry savanna, 3 % desert.
- **Essos:**
  - the Dothraki sea is grassland: 60 % temperate grassland, 20 % xeric shrubland (canon E6);
  - the Red Waste is 75 % desert (E7);
  - large rivers through dry land carry riparian woodland.
- **Haunted Forest (64–66° N):** tundra 59 %, bog 41 %.

Whole world: 102.6 Pg of NPP per year; 13.5 M km² of the land is lake water (E12: mostly an inland sea of the
south-western continent and two basins of Sothoryos). Against Earth's WWF shares the
world has more desert and xeric shrubland (33 % against 21 %) and temperate grassland (12 % against 7.5 %), and less
boreal forest and tundra (5 % and 3 % against 11 % and 9 %). That follows from its drier interiors and its smaller
northern land.

Run g is deterministic (two runs are identical) and takes ~15 s.

## Still open

- Mangroves, tropical conifer forests and the monsoon and cloud forests are not separate biomes.
- One biome per 10 km cell: no mosaic, tree cover or stand age; fire and grazing are not modelled.
- The Dornish "dry savanna" is a label of the tropical branch (Tc ≥ 12 °C). Its vegetation is a thorn scrub of
  the Sahel type.

## Sources

- Prentice, I. C., Cramer, W., Harrison, S. P., Leemans, R., Monserud, R. A., Solomon, A. M. (1992), A global biome
  model based on plant physiology and dominance, soil properties and climate, J. Biogeogr. 19, 117.
- Haxeltine, A., Prentice, I. C. (1996), BIOME3, Global Biogeochem. Cycles 10, 693.
- Kaplan, J. O. et al. (2003), Climate change and Arctic ecosystems (BIOME4), JGR 108, 8171.
- Dallmeyer, A., Claussen, M., Brovkin, V. (2019), Harmonising plant functional type distributions for evaluating
  Earth system models, Clim. Past 15, 335 — the mega-biome limits.
- Olson, D. M. et al. (2001), Terrestrial ecoregions of the world, BioScience 51, 933; areas from Bertzky, B. et
  al. (2013), Terrestrial Biodiversity and the World Heritage List, IUCN.
- Lieth, H. (1975), Modeling the primary productivity of the world (Miami model), in Primary Productivity of the
  Biosphere, Springer.
- Field, C. B. et al. (1998), Primary production of the biosphere, Science 281, 237.
- Thornthwaite, C. W., Mather, J. R. (1955), The water balance, Publ. Climatol. 8.
