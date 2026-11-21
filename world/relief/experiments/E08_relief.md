# E08 — Global relief (bedrock elevation)

## Goal

The mean elevation of every 10 km cell, land and sea floor, plus the expected relief inside the cell, built only from
the layers above it (plates, canon relief, history) so that regional maps can be derived from it. Bedrock only:
no ice, lakes or river valleys yet (those layers come later and will modify it).

## Method

**Land** = plain + ranges − canyons.

- *Plain*: a level and an amplitude for each tectonic province (shields 450 ± 250 m, platforms 250 ± 120 m,
  basins 100 ± 50 m, plateaus 1500 ± 300 m, …), blurred across province borders (σ ≈ 50 km) and varied by a
  fractal field. It falls toward the coast over a coastal plain 150 km wide on passive margins and 40 km on
  active ones.
- *Ranges*: every range has an Earth analogue whose highest summit (`peak_m`) sets the mean height of a crest cell
  (0.6 × peak above the plain). Across the range the height follows a sine profile from its edge to its crest,
  is broken into ridges and valleys by a ridged fractal field, and is blurred by σ ≈ 26 km so that the hand-drawn
  envelopes lose their corners.
  - Canon ranges (E04 envelopes, 22 features) also take the density of relief drawn on the official map
    (70–100 % of the crest height).
  - Point features (Dragonmont, Mother of Mountains, Marahai) are cones.
  - Ranges outside the canon follow the history provinces: worn old belts, collision belts, arcs, and the four
    active margins (Andes-type on the southern continent, Japan-type east of Essos and Ulthos, Sumatra-type
    on the Jade Sea).
- *Canyons*: the Great Sand Sea and the Dry Deep are cut 1.0 and 1.5 km into their plateau. The Neck is held
  near 30 m.

**Sea floor**

- *Oceanic crust*: GDH1 basement depth from E07, raised by sediment, deepened at trenches, and given abyssal
  hills (± 150 m).
  - Sediment: 5 m/Myr pelagic, at most 1 km, plus a terrigenous apron of 2 km × exp(−d/300 km) from land.
    Isostatic rise is 0.57 × thickness: (3300 − 2000)/(3300 − 1030).
  - Trenches: 3 km × exp(−(d/50 km)²).
- *Continental margin*: a shelf to 130 m over the inner 40 % of the margin, then a smooth slope down to the
  nearest oceanic floor.
- *Stepstones*: 40 m.

**Sub-grid relief** = 40 m + 0.9 × range height above the plain (+ half the canyon depth); 300 m on the sea floor.

## Runs

| Run | Change | Mean land, m | Land > 1000 / 2000 m | Mean ocean, m | Verdict |
|---|---|---|---|---|---|
| a | first version | 641 | 17 % / 1.3 % | 4515 | no: province borders were steps (mesas with straight cliffs); too little high ground; ocean too deep (no sediment) |
| b | plains blurred (seam-aware smoothing added to the grid tools); sediment; broader crests | 693 | 22 % / 2.0 % | 4228 | no: the southern continent had no Andes; E06 trenches lie far off the generated coasts, so the E07 margin belts missed the land |
| c–e | E07 run d: belts measured from the facing coast; narrow sea shoulders 200 km | 854 | 29 % / 7.4 % | 4228 | partly: Andes and arcs in place; ranges still flat-topped polygons; point volcanoes missing |
| f | point features as cones | 854 | 29 % / 7.4 % | 4228 | partly: in close-up the hill envelopes are mesas, and the Neck was 700 m high because "plateau or basin" shared one code |
| g–h | ridged noise, sine profile, blurred ranges and Neck weight | 725 | 19 % / 4.1 % | 4228 | partly: Neck still high (province code) |
| i | basins get their own province code (E07 run d) and a 100 m plain | 723 | 19 % / 4.1 % | 4228 | yes for land; three open questions (below) |
| j | summit heights verified (Etna 3326, Narodnaya 1894, Olympus 2917); Marahai as a ring island; E07 run e with Mid-Sunset at 75 km/Myr | 723 | 19 % / 4.1 % | 4189 | partly: the median now matches Earth (4701 m vs 4700 m), the mean does not |
| k | E07 run e complete: oceanic plateaus +2 km, passive margins 300 km | 723 | 19 % / 4.1 % | 4041 | partly: distribution too narrow, peak at 5.0–5.5 km |
| l | residual depth (dynamic topography), rms 500 m, zero mean | 723 | 19 % / 4.1 % | 4062 | yes then |
| m | E07 run f: new ridges (equatorial, extinct East Ulthos), Antarctic ring both sides, reach rule | 723 | 19 % / 4.1 % | 4030 | yes then; the gap is not in the ridges (below) |
| n | sediment calibrated to GlobSed; sediment thickness saved as a layer | 723 | 19 % / 4.1 % | 3953 | yes then; deep tail twice Earth's |
| o | E07 run g; sediment refitted; residual depth with and without (rms 0–500) to measure its effect | 723 | 19 % / 4.1 % | 3881–3903 | calibration run |
| p | residual depth rms 500 → 261 m | **723** | **19 % / 4.1 %** | **3880** | **yes** |

## Closing the open questions (runs j–l)

**Summit heights.** Every `peak_m` was checked against summit lists and encyclopaedias. Three were corrected:
- Etna 3357 → 3326 m (quoted values range 3326–3357 m, because the summit changes with eruptions);
- Narodnaya 1895 → 1894 m;
- Olympus 2918 → 2917 m.

The rest stand. Galdhøpiggen 2469, Ben Nevis 1345, Sgùrr Alasdair 992, Merrick 843, Brocken 1141, Signal de
Botrange 694, Grand Ballon 1424, Aneto 3404, Jabal Sawda ~3000 (official 3015, surveyed 2999), Dena 4409,
Jengish Chokusu 7439, Elbrus 5642, Puncak Jaya 4884, Emi Koussi 3415, Stromboli 924, Mount Mitchell 2037,
Ramelau 2963, Aconcagua 6961, Fuji 3776 and Kerinci 3805 are confirmed.

**Marahai.** The island was not lost: E03 kept it. It is a ring ~300 km across around a lagoon, and the E04
point sits in the lagoon, so the cone fell on water. The ring now takes a range profile on its own land (highest
cell 593 m). A 300 km ring is far larger than any caldera (Toba is ~100 km). The closest Earth form is a curved
island arc around a small sea, such as the Banda arc; Marahai lies on the E06 Jade arc.

**Ocean depth.** The comparison was made against Earth's hypsometry by 500 m bands (1° grid). The gap had three
causes:
1. The floor was too old: mean age 80 Myr against Earth's 64.2 Myr. Faster spreading of the Sunset Sea gives 67.
2. Earth's ~5 % of oceanic plateaus, which stand 2–3 km above the floor, were missing.
3. The floor had no dynamic topography. Residual depths on Earth vary by ±1 km at wavelengths of ~1000 km, and
   without them the distribution was too narrow.

Mean ocean depth went from 4228 m to 4062 m, against Earth's 3688 m.

| Depth band, m | 0–500 | 500–2500 | 2500–4000 | 4000–5000 | 5000–6000 | > 6000 |
|---|---|---|---|---|---|---|
| This world, % of ocean | 9.5 | 4.3 | 19.8 | 37.2 | 27.2 | 2.2 |
| Earth, % of ocean | 10.2 | 9.2 | 27.5 | 33.2 | 19.2 | ~0.6 |

The remaining gap is about 370 m (10 %). The shelf and the abyss now match. What is missing is the 500–4000 m
range:
- continental slopes and rises;
- ridge flanks: our ridges total ~49 000 km of axis against Earth's ~65 000 km, and three of five are slow;
- back-arc basins, island arcs and seamount chains, which the history does not contain.

**Run m: ridges added** (E07 run f). The mean depth moved only from 4062 to 4030 m. The mean age is already
Earth-like (68 against 64.2 Myr), so the remaining ~340 m does not come from the age of the floor. It comes from
what lies on the floor: Earth's continental rises carry sediment up to ~10 km thick (GlobSed), while our apron is
at most 2 km. That is the missing 500–3500 m band.

**Run n: sediment calibrated.** GlobSed (Straume et al. 2019) gives three figures:
- the continental margins, 12.9 % of the ocean, hold 3044 m of sediment on average;
- oceanic floor more than 200 km beyond them holds 404 m;
- the oceans hold 337 million km³ in all.

Sediment = min(r × age, 1 km) + A × exp(−(d/L)²), where d is the distance from land. Fitting r, A and L to the
three figures gives r = 2.77 m/Myr (within the 1–10 m/Myr of pelagic sedimentation), A = 3137 m and
L = 558 km. With these the model reproduces all three: 3044 m, 404 m and 336.9 million km³. An exponential
apron could not fit them with a positive pelagic rate.

Mean ocean depth fell to 3953 m. The sediment thickness is now a layer of its own (`sediment_m`), which later
layers can use: soils and resources, for example oil and gas in thick margin sediment.

The gap that remains is 265 m (7 %). It now lies in the deep tail: 12.5 % of our ocean is deeper than 5.5 km,
against 6.3 % on Earth. Two causes:
- the 160 Ma floor that no ridge reaches (11 % of the ocean, E07 run f);
- the negative half of the residual-depth field falling on that old floor.

Shelf, slope and abyssal bands are each within 1–3 percentage points of Earth's.

**Runs o–p: the deep tail.** Every cell deeper than 5.5 km reached that depth only through the added terms
(trenches, residual depth, abyssal hills). The age–depth curve with sediment never goes past 5.5 km. So the old
floor was not the cause; the residual-depth field was too strong.

The amplitude was set from Earth's hypsometry. With rms 261 m, 5.6 % of our ocean is deeper than 5.5 km against
6.3 % on Earth, and 0.8 % deeper than 6 km against ~0.6 %. Peaks of about ±0.8 km stay within the observed ±1 km.

The history was also brought closer to Earth (E07 run g: extinct Shivering Sea ridge, faster Antarctic ring).
The sediment was refitted to GlobSed for the younger floor: 3.25 m/Myr, 3134 m, 554 km; it again gives
3044 m / 404 m / 337.1 million km³.

Mean depth is 3880 m against 3688 m, a 192 m (5 %) gap. It now sits in the middle of the abyss: 40 % of our ocean
lies at 4.5–5.5 km against 30 % on Earth, and 3.5–4.5 km is short by 5 points. That is the age–depth curve itself
(GDH1 plus sediment for floor 40–120 Ma old) applied to our ocean, not a missing feature.

## Result of run p (land as in run i)

| | This world | Earth |
|---|---|---|
| Mean land elevation | 723 m | 840 m |
| Mean ocean depth | 3880 m | 3688 m |
| Median ocean depth | 4365 m | — |
| Ocean deeper than 5.5 km / 6 km | 5.6 % / 0.8 % | 6.3 % / ~0.6 % |
| Sediment: margins / far ocean / total | 3044 m / 404 m / 337.1 M km³ | 3044 m / 404 m / 337 M km³ |
| Shelf (sea floor above 200 m) | 26.9 M km² | 28.4 M km² |
| Highest cell mean | 5823 m (Andes-type margin) | — |
| Deepest cell | −8795 m (trench) | −10 924 m (Challenger Deep) |

- Earth's mean land elevation includes the Antarctic and Greenland ice sheets. Earth's shelf area is land plus
  shelf (177.4 M km²) minus land (149 M km²).
- The hypsometric maximum of Earth's deep-sea floor is at 4700 m.
- Coastline: every land cell is above 0 and every sea cell below 0 (the E03 outline is kept exactly).
- The Neck averages 128 m.
- Canon ranges reach 50–80 % of their analogue's summit in the highest cell mean, as expected for 10 km means
  (Bone Mountains 3464 m vs Tian Shan 7439 m; Red Mountains 2027 m vs Aneto 3404 m).

## Known limits

- The ocean is ~190 m (5 %) deeper on average than Earth's: too much floor at 4.5–5.5 km, too little at 3.5–4.5 km.
- Hill envelopes still show traces of their polygons in close-up; regional maps will need their own detail.
- Ulthos' east coast is passive: its trench lies ~2000 km offshore (as Tonga–Kermadec lies off Australia).
- Bedrock only. Ice sheets (Land of Always Winter, south polar continent), lakes and river incision come with the
  climate and hydrology layers.

## Sources

- Mean land elevation 840 m, mean ocean depth 3688 m — Britannica, "Ocean: Major subdivisions" (hypsometric curve).
- Hypsometric maxima at 100 m and 4700 m — Britannica, "Hypsometric curve".
- Land 149 M km², land plus shelf 177.4 M km² — Segar, *Introduction to Ocean Sciences*, Appendix 2 (LibreTexts).
- GDH1 — Stein & Stein (1992), Nature 359, 123; Parsons & Sclater (1977), JGR 82, 803.
- Ocean hypsometry by 500 m bands (1° grid) — S. Dutch, "Earth hypsometry" table (stevedutch.net).
- 0–200 m is 7.6 % of the ocean — Britannica, "Ocean: Major subdivisions".
- Residual depth ±1 km at ~1000 km wavelength — Hoggard, White & Al-Attar (2016), Nature Geoscience 9, 456.
- Oceanic plateaus ~5 % of the oceans, 2–3 km high — see E07 sources.
- Mean ocean floor age 64.2 Myr — Seton et al. (2020), G³ 21.
- Summit heights — Store norske leksikon (Galdhøpiggen, Aneto); Britannica (Emi Koussi, Southern Uplands);
  summit and ultra-prominence lists (ii.uib.no, Wikipedia lists of country high points and of volcanoes);
  Country Highpoints survey of Jabal Sawda (2018); Wikipedia articles for the individual peaks.
- Sediment budget (margins 12.9 % of the ocean, 3044 m; far ocean 404 m; 337 million km³) — Straume et al. (2019),
  GlobSed, G³ 20, 1756, and the University of Oslo press release (phys.org, 2019).
