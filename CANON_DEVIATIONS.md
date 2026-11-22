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

## 5. The Haunted Forest is sparse taiga at the edge of the forest-tundra

- **Canon:** the haunted forest is a vast forest beyond the Wall, from the Wall to the lands of Thenn, with ash, oak,
  chestnut, fir, ironwood, sentinel, soldier pine and weirwood; from the Wall it looks like a vast green expanse
  ([A Wiki of Ice and Fire, Haunted forest](https://awoiaf.westeros.org/index.php/The_Haunted_Forest),
  [Beyond the Wall](https://awoiaf.westeros.org/index.php/Beyond-the-Wall)).
- **Model:** at 64–66° N, just north of the Wall, the warmest month averages ~7 °C (7.6 before the ice layer, 7.0
  with it). That is below the ~10 °C that the polar tree line follows. North of the Wall the land is sparse taiga in
  sheltered valleys and lowlands, opening into forest-tundra. Dense forest, broadleaf trees included, grows only in
  the south of the North, below ~60° N: Winterfell's warmest month is ~11 °C.
- **Why:** the climate model (experiment `world/climate/experiments/E09_temperature.md`) is calibrated on 78 Earth
  stations and reproduces Scandinavian and Alaskan summers within 1–3 °C. It includes the warm drift that the owner
  chose to keep the North as mild as possible (decision 7). The cool summers come from the world's geography: the
  lands at 64–70° N lie between two oceans, with an icy polar sea and the Land of Always Winter to the north.
  Uncertainty is about ±3 °C, so even at its upper end the Haunted Forest stays at the tree line. Project owner's
  decision, 2026-11-21.

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

