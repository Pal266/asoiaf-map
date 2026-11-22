# E14 — Lithology

**Goal:** the rock at the surface of every 10 km cell, as a class of the Global Lithological Map (GLiM, level 1),
built from the geological history (E07) and from the layers above it. It is the parent material of the soils (E13)
and, later, the base of the resources.

## Research

- **GLiM** (Hartmann & Moosdorf 2012) maps Earth's land surface in 16 level-1 classes: unconsolidated, siliciclastic,
  mixed and carbonate sediments, evaporites, metamorphics, acid / intermediate / basic plutonics and volcanics,
  pyroclastics, ice and water. Land shares: about 64 % sediments (a third of them carbonates), 13 % metamorphics,
  7 % plutonics, 6 % volcanics, 10 % water and ice. Unconsolidated sediments alone cover 24.6 % of the land
  (Boerker et al. 2018).
- **Provinces and rocks** (textbook geology, e.g. Condie, *Earth as an Evolving Planetary System*):
  - shields: gneiss and granite–greenstone terranes, with small sedimentary covers and dykes;
  - platform covers on cratons: flat-lying sandstones, shales and carbonates;
  - young orogens: folded sediments (flysch, carbonate platforms) with metamorphic cores and granites;
  - arcs: andesites, pyroclastics and granodiorite batholiths;
  - rifts and passive margins: rift sediments and flood basalts;
  - large igneous provinces: flood basalts;
  - intracontinental basins: thick young sediments, evaporites in closed basins.
- **Surface deposits:**
  - alluvium along large rivers and in deltas and marshes;
  - evaporites in salt pans;
  - aeolian sand seas on low desert plains (about a fifth to a third of the hot deserts, Pye & Tsoar 2009);
  - loess on semi-arid to subhumid plains downwind of deserts and outwash plains, typically up to several hundred km
    from the source (Pye 1995).

## Method

1. **Bedrock by province.** Each province of E07 has an assumed composition (`province_composition` in
   `data/lithology.json`). One smooth uniform field (a ~130 km octave plus a ~30 km octave, rank-mapped to a uniform
   distribution) picks the class of each cell from the cumulative shares, so rocks form coherent outcrops in the
   right proportions.
2. **Shield and cover.** Cratons crop out as shield where local relief is above 250 m or where a second field says so
   (at least 20 % of a craton); the rest is platform cover.
3. **Warped contacts.** The history's provinces are polygons and buffers. The province is read at a point displaced
   by a smooth random warp of ~0.8°, so contacts wander like real ones.
4. **Canon volcanic relief:** intermediate volcanics and pyroclastics.
5. **Deposits:**
   - alluvium: discharge ≥ 100 m³/s on relief below 150 m, and marshes;
   - evaporite: salt flats;
   - aeolian sand: warm (MAT ≥ 5 °C) arid plains (P < 250 mm, relief < 150 m), on 35 % of them;
   - loess: low plains (relief < 200 m, coldest month < 5 °C, P < 800 mm, not arid) within ~800 km of desert dust.
6. Lakes become water and ice sheets become ice.

## Experiments

| Run | Change | Sediments | Metam. | Plut. | Volc. | Unconsol. | Note |
|---|---|---|---|---|---|---|---|
| b | compositions with the shield / cover split | 0.704 | 0.132 | 0.080 | 0.084 | 0.229 | polar "sand seas" on cold arid plains |
| c | aeolian sand and dust need MAT ≥ 5 °C | 0.700 | 0.132 | 0.082 | 0.087 | 0.223 | basalt bands along the rifted margins |
| d | warped contacts; rifted margins mostly sedimentary | 0.705 | 0.135 | 0.083 | 0.077 | 0.224 | final |
| GLiM | | 0.71 | 0.144 | 0.078 | 0.067 | 0.27 | shares of land without water and ice |

Carbonates make up 26 % of the sediments (GLiM about a third). Deposits (share of all cells' land): alluvium 3.7 %,
evaporite 1.5 %, aeolian sand 7.4 %, loess 8 %.

Run d is deterministic (two runs are identical).

## World (run d)

- **Westeros south of the Wall:**
  - shares: deposits 24 %, siliciclastic 22 %, carbonate 22 %, mixed 15 %, metamorphic 8 %, granites 5 %;
  - the North: shield gneiss and granites around Winterfell and in the northern mountains, alluvium in the Neck;
  - the Reach and the Stormlands: platform sandstones and limestones, with loess on the dry plains;
  - Dorne: aeolian sand and limestone, with a metamorphic and granite core in the Red Mountains.
- **Beyond the Wall:** platform sediments with basalts on the old margins.
- **Western Essos:**
  - deposits 37 %: the loess belt of the Dothraki Sea and sand in the deserts;
  - Valyria and the Fourteen Flames: pyroclastics and andesites.

## Soils (E13 with lithology)

The soil key now takes parent material from this layer (see E13, addendum):

| Group (M ha) | without lithology | with lithology |
|---|---|---|
| Arenosols | 1470 | 895 |
| Vertisols | 294 | 74 |
| Andosols | 33 | 176 |
| Nitisols | 0 | 21 |
| Leptosols | 1380 | 1518 |
| Luvisols | 609 | 688 |

The Earth test is unchanged (52 % of sites, area log RMSE 0.26): Earth has no lithology input, and the key's old
rules apply when the lithology is absent.

## Still open

- The compositions are assumed per province, not derived from a simulated rock history (sedimentation, burial,
  exhumation).
- No age of rocks, no faults, no ore bodies; the resource layer will need them.
- Glacial till and outwash are not mapped: the ice layer gives only today's sheets.

## Addendum (2026-11-22, from E17): terranes and point volcanoes

- **Terranes with a composition of their own** (`terranes` in `lithology.json`). A terrane is a history block
  (tested in true coordinates) or a canon feature. It is drawn with the ~30 km field only, in its own class
  sequence. The fixed class order of the provinces lets only neighbours in the order touch, so limestone never meets
  granite; a terrane's sequence can put them side by side.
  - Iron Islands (decision 12): granite and diorite bosses in slate, limestone and greenstone (Cornwall).
  - Westerlands hills (decision 13): turbidites, metavolcanics, greenstones, felsic volcanics and granitoids
    (Lachlan orogen).
- **Point volcanoes.** Volcanoes drawn as points on the map (Dragonmont) are not in the feature raster. The land
  within 30 km of the point is now volcanic, so Dragonstone has andesite and obsidian.
- **Effect.** The world's lithology shares are unchanged to three decimals.

## Sources

- Hartmann, J., Moosdorf, N. (2012), The new global lithological map database GLiM, G3 13, Q12004.
- Boerker, J., Hartmann, J., Moosdorf, N., Amann, T. (2018), Global distribution of unconsolidated sediments,
  Geoderma 313, 152.
- Condie, K. C. (2016), Earth as an Evolving Planetary System, 3rd ed., Academic Press.
- Pye, K. (1995), The nature, origin and accumulation of loess, Quat. Sci. Rev. 14, 653.
- Pye, K., Tsoar, H. (2009), Aeolian Sand and Sand Dunes, Springer.
