# Canon deviations

## 1. Regular one-year seasons

- **Canon:** seasons last for years, sometimes decades, and their length cannot be predicted. The year itself is as
  long as Earth's ([A Wiki of Ice and Fire, Measurement](https://awoiaf.westeros.org/index.php/Month)).
- **Model:** one year of four regular seasons, as on Earth.
- **Why:** the planet has Earth's orbit and axial tilt (`world/planet/PLANET.md`). With these, seasons follow the
  tilt of the axis over one orbit, and last one year. Canon gives no natural cause for its long, unpredictable
  seasons, and every proposed explanation needs a different planet or star system (an unstable orbit or axis, a second
  star). Rule 1 requires the scientifically grounded option.

## 2. Westeros is not the size of South America

- **Canon:** George R. R. Martin has said that Westeros is similar in size to South America; he has also said that the
  Wall is 100 leagues (300 miles) long ([A Wiki of Ice and Fire, Westeros](https://awoiaf.westeros.org/index.php/Westeros)).
  The two statements disagree.
- **Model:** Westeros is scaled by the Wall (~480 km), which gives ~4 800 km north to south and ~1 450 km at its
  widest. South America is ~7 500 km long and ~17.8 million km² in area.
- **Why:** the Wall's length is a stated measurement that fixes the scale of the canon maps; the comparison with South
  America is a loose remark that no map supports. Project owner's decision (`world/canon/CONSTRAINTS.md`, section 7).

## 3. Directions on the map are not true in eastern Essos and beyond the Wall

- **Canon:** the map of the known world is drawn with north up everywhere, and the books' directions ("north of",
  "east of") follow it.
- **Model:** the map lies on the globe as a sinusoidal projection centred on the Wall, so local north leans away from
  the map's vertical with distance from the Wall's meridian. South of the Wall, Westeros keeps its directions within
  11° (95 % of the land); beyond the Wall within 31°; in western Essos within 25°; in central Essos within 42°; east of
  the Bone Mountains within 51°. Distances along the parallels and all areas stay as drawn.
- **Why:** a flat map 11 800 × 7 800 km cannot keep scale, shape and direction at once on a sphere. Of the readings
  tested, only the sinusoidal one keeps the Wall's length, the outlines of Westeros, all areas and the canon climates
  together; it gives up directions where canon describes the least. Project owner's decision
  (`world/canon_outline/experiments/E02_projection.md`).

## 4. Essos is the size of Asia, not of Eurasia

- **Canon:** Essos is "roughly the size of Eurasia" (~55 million km²), a remark outside the books (grade B)
  ([A Wiki of Ice and Fire, Across the narrow sea](https://awoiaf.westeros.org/index.php/Across_the_narrow_sea)).
- **Model:** Essos covers ~41 million km² (~21 on the map, ~20 east of it), about the size of Asia (44.6 million km²).
- **Why:** on the globe the map's east edge lies at 90–147° E across the latitudes of Essos (15–55° N), and its east coast
  can lie no farther than ~165–168° E if the Sunset Sea is ~180° wide (`world/canon/CONSTRAINTS.md`, decision 5). That
  leaves room for ~24 million km² east of the map. Project owner's decision (decision 6, revised 2026-11-18).

## 5. The Haunted Forest is tundra; the North is forest-tundra and bog, with taiga in its south

- **Canon:** the haunted forest is a vast forest beyond the Wall, from the Wall to the lands of Thenn, with ash, oak,
  chestnut, fir, ironwood, sentinel, soldier pine and weirwood; from the Wall it looks like a vast green expanse
  ([A Wiki of Ice and Fire, Haunted forest](https://awoiaf.westeros.org/index.php/The_Haunted_Forest),
  [Beyond the Wall](https://awoiaf.westeros.org/index.php/Beyond-the-Wall)). The North has large forests, such as
  the wolfswood around Winterfell.
- **Model:**
  - at 64–66° N, just north of the Wall, the warmest month averages ~7 °C and summers have ~90 degree-days above
    5 °C, a quarter of what trees need. The vegetation layer (E15) makes the Haunted Forest tundra (59 %) and bog
    (41 %); trees grow there only as scattered krummholz in sheltered valleys;
  - south of the Wall the North is open:
    - 60–63.75° N: tundra 73 %, bog 22 %;
    - 56–60° N: bog 35 %, tundra 31 %, forest-tundra 30 %;
    - 52–56° N: bog 39 %, forest-tundra 38 %, taiga 15 %;
  - Winterfell (warmest month ~11 °C) stands in bog and forest-tundra, the open, park-like edge of the taiga;
  - closed forest begins at the Neck.
- **Why:**
  - the climate model (experiment `world/climate/experiments/E09_temperature.md`) is calibrated on 78 Earth
    stations. It includes the warm drift that the owner chose to keep the North as mild as possible (decision 7);
  - re-examined for the vegetation layer (E09, addendum 2): at Earth's maritime subarctic stations the model's
    summers have no cold error (mean +0.1 °C north of 55° N). The cool summers come from the world's geography.
    The seas on both sides of the North are cold (July at 56° N: 10 °C to the west, 5 °C to the east), and the polar
    sea keeps ice into July. The land does not warm above the sea. On Earth, forest at 64–66° N grows only in
    continental interiors with warm summers, not on maritime coasts such as Iceland's;
  - the forest limit (350 degree-days above 5 °C) is the BIOME1 value (Prentice et al. 1992).
- **Decision:** project owner, 2026-11-21 (sparse taiga, after E09), extended 2026-11-22 after E15. The owner
  preferred this to moving the map south.

## 6. The Land of Always Winter has no ice sheet

- **Canon:** the Land of Always Winter is the northernmost part of Westeros, permanently frozen (grade A)
  ([A Wiki of Ice and Fire, Land of Always Winter](https://awoiaf.westeros.org/index.php/Lands_of_Always_Winter)).
- **Model:** it is polar tundra on permafrost, with a small ice cap (~28 000 km²) in its mountains at ~72° N. Its
  warmest month averages ~2.7 °C; snow lies for most of the year. "Permanently frozen" holds for the ground
  (permafrost), not as an ice sheet.
- **Why:** in the ice model (experiment `world/ice/experiments/E11_ice.md`) the northern ice is bistable. Either the
  sheet spreads over the whole land beyond the Wall, into the mountains south of it, or it does not form at all. No
  parameter value confines it to the Land of Always Winter. That is the small ice-cap instability of energy-balance
  models. Of the two states the owner chose the one that keeps the haunted forest and the North free of ice. Its
  mass-balance parameters are within the literature range, but on Earth they give ~34 % too little ice. Project
  owner's decision, 2026-11-22.

## 7. The Reach and the Crownlands have a dry Mediterranean climate; their rivers are small

- **Canon:** the Reach is the most fertile region of Westeros; the Mander is one of its major rivers, and the
  Blackwater Rush reaches the sea at King's Landing (the official map; grade A)
  ([A Wiki of Ice and Fire, The Reach](https://awoiaf.westeros.org/index.php/Reach),
  [Mander](https://awoiaf.westeros.org/index.php/Mander)).
- **Model:** the south of Westeros gets ~250–350 mm of precipitation a year, falling in winter, with dry summers.
  King's Landing has 268 mm and Oldtown 325 mm. The canon rivers flow along their drawn courses (decision 9), but
  the Mander carries ~240 m³/s from 600 000 km² and the Blackwater Rush ~100 m³/s. The Honeywine and the Greenblood
  are seasonal streams. The Reach's fertility is that of California's Central Valley or Andalusia: rich soils and
  river water, and the vineyards of the Arbor.
- **Why:** a west coast at 30–40° N on a Pacific-sized ocean (decision 5) lies under the subtropical high, as do
  California, central Chile and Morocco. Recalibrating the precipitation model with Earth's river discharge did not
  change it (King's Landing 268 → 269 mm; experiment `world/hydrology/experiments/E12_hydrology.md`). Project
  owner's decision, 2026-11-22.

