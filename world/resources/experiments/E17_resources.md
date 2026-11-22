# E17 — Ore bodies and mineral resources

**Goal:** the mineral deposits of the world (their type, place, size, grade and contained metal) and the area
resources a medieval economy uses: coal, salt, bog iron, peat, clay, building stone, obsidian, sulfur. Built from
the lithology (E14), the geology (E16: ages and faults), the climate, the soils (E13) and the vegetation (E15).

## Research

- **Three-part assessment** (Singer & Menzie 2010). USGS estimates undiscovered resources in three steps:
  - delineate the tracts permissive for each deposit type of the descriptive models (Cox & Singer 1986);
  - estimate the number of deposits from deposit densities;
  - take their sizes and grades from grade-and-tonnage models.

  Run forward, the same three steps make a world.
- **Deposit density** (Singer 2008, general across 10 deposit types and 108 control tracts):
  log10(deposits per 100,000 km²) = 4.2474 − 0.5146 log10(tract area, km²) − 0.2254 log10(median tonnage, Mt).
  For porphyry copper, densities range from 2 to 128 deposits per 100,000 km², with a median of 15.
- **Grade-and-tonnage models:**
  - median tonnage and its spread (sd of log10 tonnage) for all models: USGS OFR 90-412, Table 1;
  - grades of porphyry copper and the skarn and replacement models: OFR 2007-1214, Table 5.1-1;
  - grades of felsic and mafic massive sulfides: OFR 2009-1034, Table 2;
  - sediment-hosted copper (reduced facies): Cox et al. (2003), via OFR 2007-1214;
  - low-sulfide gold-quartz veins (Phanerozoic and Archean tonnages, a gold grade of ~14 g/t, Au:Ag 5:1):
    OFR 2002-195K;
  - no grade was found for epithermal veins, placers, MVT, SEDEX, tin, iron formations, ironstones, lateritic
    nickel and mercury. These grades are typical values and are marked as approximate in `resources.json`.
- **Geological settings** (Cox & Singer 1986; Leach et al. 2005 for sediment-hosted lead-zinc):
  - porphyries and epithermal veins form in arcs;
  - skarns form at intrusive contacts with carbonates;
  - orogenic gold forms on faults in metamorphic belts;
  - MVT deposits form in foreland carbonate platforms, from brines driven by orogenic topography;
  - SEDEX deposits form in rift basins, at their bounding faults;
  - tin forms in the granites of collisional belts;
  - banded iron formations are Archean to Paleoproterozoic;
  - oolitic ironstones are Phanerozoic marine sediments;
  - nickel laterites form on mafic rock under a humid tropical climate.
- **Earth references:**
  - USGS Mineral Commodity Summaries 2026: copper identified 2.1 Gt (past production included) plus 3.5 Gt
    undiscovered; zinc identified ~1.9 Gt; lead > 2 Gt; nickel > 350 Mt; mercury ~600 000 t;
  - well-explored deposits in the USGS world databases: 422 porphyry copper (OFR 2008-1155) and 867 massive
    sulfide deposits (OFR 2009-1034).

## Method (scripts/resources.py)

1. **Tracts.** Each of 21 deposit models has a tract rule (`resources.json`). A rule combines:
   - tectonic province, surface rock or bedrock, and age;
   - nearness to a rock (intrusions, evaporites), to faults of given types, or to a province (forelands);
   - climate (laterites);
   - for placers, alluvium within 150 km of a lode-gold tract.

   Outcrops closer than 30 km form one tract, as USGS tracts are drawn at regional map scale.
2. **Number of deposits.** Each tract of at least 200 km² gets Singer's median number for its area and the model's
   median tonnage; the actual number is a Poisson draw. Deposits go to random cells of the tract, and vein types
   prefer cells within 20 km of a fault.
3. **Size.** Tonnage is lognormal, from the model's median and spread, truncated at ±3 standard deviations. Grades
   are the model's medians, and contained metal is tonnage × grade.
4. **Area resources:**
   - coal measures: sandstone and shale of the two great coal ages, Carboniferous–Permian and Cretaceous–Paleogene;
   - salt: evaporites and salt flats;
   - bog iron: Histosols and Gleysols outside the tundra;
   - peat: Histosols;
   - clay: alluvium, shale and Vertisols;
   - granite, limestone, marble (limestone at intrusive contacts), sandstone, slate;
   - obsidian: volcanics younger than 5 Ma;
   - sulfur: volcanics younger than 1 Ma.

## Runs

| Run | Change | Deposits | Zn, Gt | Sn, Mt | Note |
|---|---|---|---|---|---|
| a | first rules | 69 070 | 38.7 | 762 | permissive outcrops counted as separate small tracts |
| b | outcrops closer than 30 km form one tract | 50 324 | 31.5 | 563 | tracts now of USGS size (~10⁴–10⁵ km²); the excess comes from too much permissive area |
| c | MVT in forelands (≤ 300 km from an orogen), SEDEX at basin faults, tin within 100 km of thrusts and sutures | 43 823 | 16.4 | 322 | |
| d | full chain rerun from the history: the Iron Islands terrane (decision 12), volcanic Dragonstone | 43 675 | 17.0 | 325 | Iron Islands: tin and some gold, but no iron or lead |
| e | Iron Islands: own composition and outcrop sequence; polymetallic veins added | 68 724 | 21.3 | 335 | iron, lead and tin on the islands |
| f | Westerlands gold belt (decision 13); tonnage truncated at ±3σ | **69 525** | **18.5** | **289** | final |

- **Iron Islands, runs d and e.** The lithology's outcrop pattern has a scale of ~130 km, and it lays out the rock
  classes along one random field in a fixed order. So limestone never touches granite, and on small islands they
  fall on different islands: no skarns, no iron or lead.

  The islands' terrane now has its own composition: granite and diorite bosses in slate, limestone and greenstone,
  as in Cornwall. It is drawn with the ~30 km field and a sequence that puts the intrusions next to the limestone.
  Iron skarns need intermediate intrusions (Cox & Singer), hence the diorites. Polymetallic veins in any rock
  around granites were added: the main medieval source of lead and silver.
- **Westerlands, run f.** Before this, the Westerlands had 5.7 % of Westeros's gold on 9.3 % of its area: every old
  fold belt was treated alike, and orogenic gold lay only along the sparse province contacts. The Westerlands
  hills are now an accreted Variscan island-arc terrane (turbidites, greenstones, metavolcanics, felsic volcanics,
  granitoids; like the Lachlan orogen of Victoria). The geology adds a dense network of shear zones (6 300 km of
  inactive thrusts), so the whole terrane is permissive for orogenic gold and massive sulfides.
- **Tonnage truncation, run f.** Without it, one placer of 120 000 Mt (36 000 t of gold) was drawn in the Riverlands.

## Earth checks

| Quantity | World | Earth |
|---|---|---|
| Porphyry copper deposits | 878 | 422 well-explored (plus undiscovered) |
| Massive sulfide deposits | 2 770 | 867 well-explored |
| Copper | 2.66 Gt | 2.1 Gt identified + 3.5 Gt undiscovered |
| Zinc | 18.5 Gt | ~1.9 Gt identified |
| Lead | 7.9 Gt | > 2 Gt identified |
| Nickel | 0.32 Gt | > 0.35 Gt |
| Mercury | 1.65 Mt | ~0.6 Mt |

The world's land (147 M km²) is the size of Earth's.
- Copper, nickel and the deposit counts are within a factor of 2–3. Counts are expected to exceed the
  well-explored ones, because undiscovered deposits are of the same order.
- Lead and mercury are 3–4 times Earth's identified resources.
- **Zinc is about 10 times too high.** Three-quarters of it comes from MVT deposits, whose tonnage model (median
  35 Mt) was read from an OCR'd table and whose grade is approximate.
- Gold (395 000 t), silver (7.1 Mt), iron (1.24 × 10¹² t) and tin (289 Mt) have no Earth figure here. Tin is
  probably too high: Earth's tin reserves are only > 6 Mt.

## World (run f)

| Model | Deposits | Ore, Mt | Contained metal, t |
|---|---|---|---|
| porphyry copper | 878 | 444 180 | Cu 1.95 × 10⁹, Ag 5.3 × 10⁵, Au 3.4 × 10⁴ |
| epithermal gold-silver | 2 432 | 10 771 | Au 8.1 × 10⁴, Ag 1.2 × 10⁶ |
| copper skarn | 1 602 | 9 188 | Cu 1.6 × 10⁸ |
| iron skarn | 907 | 91 672 | Fe 4.6 × 10¹⁰ |
| zinc-lead skarn | 50 | 266 | Zn 1.6 × 10⁷, Pb 7.4 × 10⁶ |
| polymetallic replacement | 2 627 | 16 678 | Pb 8.7 × 10⁸, Zn 6.5 × 10⁸, Ag 2.5 × 10⁶ |
| polymetallic veins | 25 361 | 3 739 | Pb 1.9 × 10⁸, Zn 1.5 × 10⁸, Ag 5.6 × 10⁵ |
| felsic massive sulfide | 1 619 | 25 627 | Cu 3.1 × 10⁸, Zn 8.2 × 10⁸, Ag 6.4 × 10⁵ |
| mafic massive sulfide | 1 151 | 4 549 | Cu 7.7 × 10⁷ |
| orogenic gold | 8 895 | 6 570 | Au 9.2 × 10⁴ |
| Archean lode gold | 1 025 | 6 664 | Au 9.3 × 10⁴ |
| placer gold | 3 246 | 268 563 | Au 8.1 × 10⁴ |
| MVT zinc-lead | 1 167 | 347 638 | Zn 1.4 × 10¹⁰, Pb 5.2 × 10⁹ |
| SEDEX zinc-lead | 808 | 52 933 | Zn 3.0 × 10⁹, Pb 1.5 × 10⁹ |
| sediment-hosted copper | 51 | 7 320 | Cu 1.7 × 10⁸ |
| tin greisen | 1 805 | 52 992 | Sn 1.6 × 10⁸ |
| tin veins | 5 170 | 12 968 | Sn 1.3 × 10⁸ |
| banded iron formation | 1 166 | 1 550 241 | Fe 8.1 × 10¹¹ |
| oolitic ironstone | 2 852 | 1 116 355 | Fe 3.9 × 10¹¹ |
| lateritic nickel | 210 | 22 893 | Ni 3.2 × 10⁸ |
| hot-spring mercury | 6 503 | 551 | Hg 1.65 × 10⁶ |

Area resources (M km²): coal 8.8, salt 2.3, bog iron 3.9, peat 1.2, clay 24.5, granite 12.8, limestone 28.0,
marble 5.5, sandstone 22.4, slate 23.4, obsidian 6.8, sulfur 2.6.

**Canon** (CONSTRAINTS WE11, WE12, E9):
- **The Westerlands** (the hills and 100 km around): gold 3 090 t, 23 % of Westeros's on 9.3 % of its area; silver
  142 000 t, 56 % of Westeros's.
- **The Iron Islands:** iron 1.35 × 10⁸ t (7 iron skarns, ironstone), lead 1.35 × 10⁷ t, tin 1.2 × 10⁶ t, copper
  2 × 10⁶ t; gold only 89 t.
- **The Fourteen Flames:** gold 9 470 t, silver 214 000 t, copper 8.1 × 10⁷ t, tin 3.1 × 10⁵ t, iron 8.2 × 10⁹ t,
  as E9 lists them; also mercury.
- **Dragonstone:** obsidian on 29 % of the land within 30 km of Dragonmont, and two small epithermal gold-silver
  vein fields.

**Regions:**

| Region | Gold, t | Silver, t | Copper, Mt | Iron, Mt | Tin, Mt | Lead, Mt |
|---|---|---|---|---|---|---|
| Westeros south of the Wall | 13 600 | 254 000 | 196 | 24 300 | 14 | 644 |
| Beyond the Wall | 2 660 | 36 000 | 1.3 | 1 970 | 1.3 | 254 |
| Western Essos | 26 300 | 552 000 | 215 | 25 300 | 18 | 360 |

The North has wide bog iron and coal measures; the south of Westeros has tin, lead and gold belts; southern Essos
has porphyry copper, epithermal gold and silver, and mercury in its arcs.

Run f is deterministic (two runs give identical deposits) and takes ~20 s.

## Upstream changes

The two canon decisions changed the history (one block, 184 island cells) and the lithology (the two terranes and
the point volcano). The geology gained the Westerlands' shear zones. The whole chain was rerun:

| Layer | Change |
|---|---|
| Relief | ≤ 65 m on the islands |
| Temperature | ≤ 0.3 °C |
| Precipitation | ≤ 3 mm |
| Lakes | none |
| Soils of Westeros | ≤ 0.1 percentage points |
| Vegetation, lithology shares, hazard | unchanged |

## Still open

- **Zinc and tin are too high.** See the Earth checks.
- **Grades of nine models are approximate.** They are not from cited tables.
- **Regional totals vary from run to run.** Totals over a few hundred thousand km² rest on a few large deposits.
  Before the tonnage truncation the Westerlands' gold went from 430 to 3 200 t with the random sequence alone.
- **Contacts depend on the class order.** The lithology lays out rock classes along one field in a fixed order, so
  contacts between classes far apart in the order (limestone and granite) are rare everywhere except in the
  terranes with their own sequence. This undercounts skarns.
- **No fuels or other industrial minerals.** No oil, gas, phosphate, gemstones or building sand beyond the area
  resources.

## Sources

- Cox, D. P., Singer, D. A. (eds.) (1986), Mineral deposit models, USGS Bulletin 1693.
- Singer, D. A. (2008), Mineral deposit densities for estimating mineral resources, Math. Geosci. 40, 33.
- Singer, D. A., Menzie, W. D. (2010), Quantitative Mineral Resource Assessments: An Integrated Approach, Oxford
  University Press.
- USGS Open-File Report 90-412, Application of grade and tonnage deposit models (Table 1).
- USGS Open-File Report 2007-1214, chapters 5.0 and 7.0 (median tonnages and grades).
- Mosier, D. L., Berger, V. I., Singer, D. A. (2009), Volcanogenic massive sulfide deposits of the world, USGS
  OFR 2009-1034.
- USGS Open-File Report 2002-195, chapter K (low-sulfide gold-quartz veins).
- Singer, D. A., Berger, V. I., Moring, B. C. (2008), Porphyry copper deposits of the world, USGS OFR 2008-1155.
- Leach, D. L. et al. (2005), Sediment-hosted lead-zinc deposits: a global perspective, Economic Geology 100th
  Anniversary Volume.
- USGS Mineral Commodity Summaries 2026 (copper, zinc, lead, nickel, mercury, tin).
- A Wiki of Ice and Fire: Westerlands, Iron Islands, Fourteen Flames (canon).
