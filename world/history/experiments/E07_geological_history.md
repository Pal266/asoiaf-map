# E07 — Geological history and its present-day imprint

## Goal

Give every cell of the global grid the "cause" data that later layers (relief, soils, resources) need:
crust type, tectonic province, age of the last tectonic event, and the age and basement depth of the ocean floor.
The history has to explain the canon features of E04/E05 with Earth analogues and end in the E06 plates.

## Method chosen

A full backward plate reconstruction over 1 Gyr is not justified at this stage: the unknown lands are only
frames, so a kinematic model would add detail that nothing constrains. Instead E07 writes a **dated
scenario** (data/history.json, `timeline`) and paints its **present-day imprint**:

1. Ocean floor age: each ocean cell takes the distance to the nearest spreading ridge that bounds its own
   plate, divided by that ridge's half-rate, capped at the age the ridge opened.
2. Depth: basement depth from age by the GDH1 plate model (Stein and Stein 1992).
3. Continental crust: land plus a continental margin (250 km passive, 60 km near trenches).
4. Provinces: geological blocks (cratons, platforms, old orogens) with outlines bent by a smooth fractal warp;
   the unknown continents get craton cores welded by 600 Ma mobile belts (45% of the area, as Africa's
   Pan-African belts); belts along plate margins (arcs, rifted margins, collision zones); and the canon
   features with their E05 interpretation, which always win.

## Scenario (summary; ages in Ma)

3000 oldest cratons — 1800 northern shield and Dothraki craton — 600 an older supercontinent breaks up,
proto-narrow sea opens (Iapetus) — 420 northern closure, Caledonian-type belt (northern mountains,
Frostfangs core, Skagos, Ibben) — 300 last supercontinent assembled: Variscan-type Mountains of the Moon
and western Essos hills, Bone Mountains suture (Urals) — 180 breakup, Sunset Sea opens (central Atlantic) —
130 southern continent leaves Sothoryos — 100 Ulthos starts north (India/Australia) — 60 Dorne against
Westeros, Red Mountains (Pyrenees) — 30 southern Essos collision belt and plateau uplift — 25 narrow sea
reopens on the old suture (Red Sea, Wilson cycle) — 20 Bone Mountains lifted again (Tian Shan) — 10 Bleeding
Sea rift (East African Rift) — 5 Ulthos meets the Jade arcs (Banda arc), Valyrian arc — 0.4 the Doom (Santorini).

Spreading half-rates: narrow sea 10 km/Myr, Mid-Sunset 75 (full rate 150 mm/yr, in Earth's "fast" class of
75–180 mm/yr), South Sunset 15, Sothoryos–Ulthos 25, Antarctic ring 20.

## Runs

| Run | Change | Continental crust | Mean / max floor age, Ma | Mean basement depth, m | Verdict |
|---|---|---|---|---|---|
| a | first version, Parsons–Sclater depth | 34.4% | 84 / 180 | 5153 | no: each cell took the **youngest** capped age of any ridge on its plate, so distant floor got the narrow sea's 25 Ma; straight block outlines; unknown continents one uniform craton |
| b | nearest ridge per plate; fractal warp; craton cores and belts; shelves 200/60 km | 37.2% | 90 / 180 | 5259 | partly: shapes right; floor too old and deep because the Mid-Sunset ridge at 40 km/Myr needs 250 Myr to reach the far trenches (capped at 180 over a large area) |
| c | Mid-Sunset 55 km/Myr; GDH1 depth; passive margin 250 km; northern shield widened west | 38.7% | 80 / 180 | 5023 | yes for the ocean; margins found wrong later (run d) |
| d | (found while building relief, E08) margin belts measured from the coast facing the boundary, not from the boundary line; narrow sea rift shoulders 450 → 200 km; code 7 split into plateau (7) and basin (8) | 38.7% | 80 / 180 | 5023 | yes then; superseded by run e |
| e | (closing the E08 ocean-depth gap) Mid-Sunset 55 → 75 km/Myr; passive margin 250 → 300 km; oceanic plateaus (5.11 % of the ocean) | 40.1% | 67 / 149 | 4943 | yes then; superseded by run f |
| f | (owner's request: add ridges) Equatorial ridge, extinct East Ulthos back-arc ridge, Antarctic ring dates both sides, South Sunset 15 → 25 km/Myr; reach rule | 40.1% | 68 / 160 | 4904 | yes then; superseded by run g |
| g | (owner: fix the deep tail) extinct Shivering Sea ridge; Antarctic ring 20 → 32 km/Myr; floor counted by spreading class | **40.1%** | **60 / 160** | **4834** | **yes** |

Earth reference: continental crust ~41% of the surface, oldest ocean floor ~180 Ma, typical ocean floor
~4.9 km below sea level. Runs c and d are within a few percent; the remaining gap in depth is expected because
sediment and seamounts are not yet added. All 30 mapped canon features have the province E05 gives them.
Narrow sea floor is at most 25 Ma old.

Run d fixes what E08 exposed: the schematic E06 trenches lie 100–900 km off the generated coasts, so a belt
"700 km from the trench" barely touched the land and the southern continent had no Andes. Belts now start at
the coast that faces the boundary (coast cells within 1000 km of it). Active margins grew from 6.2 to 26.3
million km², rift margins are 7.5 million km². The Neck and the Shrinking Sea were coded "plateau or basin",
which gave the Neck a 1500 m plain in E08; basins now have their own code.

Run e answers the ocean-depth gap found in E08. Earth's mean ocean floor age is 64.2 Myr (Seton et al. 2020),
run d gave 80 Myr: the Pacific-sized Sunset Sea spread too slowly. At 75 km/Myr the mean is 67 Myr; the oldest floor
(149 Ma) is younger than Earth's 180 Ma because more of the Sunset floor has been subducted, which the 180 Ma
opening date allows. Oceanic plateaus (province 9) cover 5.11 % of the ocean, as on Earth (184 plateaus,
~18.5 million km²); they are placed by a seeded fractal field on floor older than 30 Ma, more than 300 km from
trenches and land. Wider passive margins bring continental crust to 40.1 % (Earth ~41 %).

Run f adds ridges at the owner's request and fixes how ridges date the floor:
- **Equatorial ridge** (active since 40 Ma, half-rate 30 km/Myr), between the Westeros and southern-continent
  plates; the analogue is the Galápagos spreading centre. E06 changes the boundary type only.
- **East Ulthos back-arc ridge**, extinct: it spread 100–45 Ma at 25 km/Myr inside the Ulthos plate, between
  Ulthos and the East Ulthos trench. Analogue: the Tasman Sea, which opened 84–52 Ma and then stopped.
- **Antarctic ring.** It now dates the floor on both sides, since every southern plate borders it. Before, it
  dated only the polar plate.
- **South Sunset ridge** at 25 km/Myr, still in Earth's slow class.
- **Reach rule.** A ridge dates only floor within half-rate × opening age of its axis. Before, every cell took
  its nearest ridge, so the narrow sea ridge (25 Ma) dated the ocean west of Westeros and the whole Shivering Sea
  as young floor. Floor no ridge reaches predates them all and is set to 160 Ma: the Shivering Sea (as the
  Amerasian basin of the Arctic) and the Summer Sea (as the eastern Mediterranean). Together that is 11 % of the
  ocean.

The mean age is 68 Myr (Earth 64.2). The new ridges make young floor, and the reach rule makes some older floor;
the mean barely moves. Its distribution is now Earth-like in where young and old floor lie.

Run g. Earth's preserved ocean floor by spreading class (Seton et al. 2020, full rates) is ultraslow 5 %, slow 33 %,
intermediate 20 %, fast 39 % and super-fast 3 %. Run f had slow 33 %, intermediate 6 %, fast 50 %, and 11 % of floor
older than every ridge (no Earth counterpart). Two changes:
- an **extinct Shivering Sea ridge** (90–35 Ma, half-rate 30 km/Myr), as the Labrador Sea and the Canada Basin
  were opened by ridges that later died;
- the **Antarctic ring** at full rate 64 mm/yr, the intermediate class of the Southeast Indian Ridge.

Now the floor is slow 20 %, intermediate 27 %, fast 47 %, with 6 % older than every ridge. That remainder is the
Summer Sea, kept old as the eastern Mediterranean is. The mean age is 60 Myr (Earth 64.2).

## Known limits

- This is a scripted history, not a simulated one: no plate motion is integrated, so the ages of continental
  provinces are assigned by block, not derived.
- Age jumps across plate boundaries are sharp; on transforms that is real (fracture zones), elsewhere it is
  an artefact of the schematic E06 boundaries.
- GDH1 coefficients are the commonly quoted fitting equations; the paper itself was not opened in this session.
- Hot-spot tracks and submerged microcontinents are not modelled; oceanic plateaus are placed statistically, not by plume history.

## Conclusion

Run g is accepted as the history layer. It feeds the relief layer next: crust type and floor age give the
base elevation, province and event age give how high and rugged the land can be.

## Sources

- Parsons, B. & Sclater, J. G. (1977) JGR 82, 803 — age–depth relation.
- Stein, C. A. & Stein, S. (1992) Nature 359, 123 — GDH1 plate model.
- Oldest ocean floor ~180 Ma; spreading rates 2–15 cm/yr — physical geology textbooks (OpenText BC, Physical Geology 2e).
- Continental crust ~41% of Earth's surface — Wikipedia, "Continental crust".
- Galápagos spreading centre and Tasman Sea opening (84–52 Ma) — standard marine geology; Gakkel ridge and Amerasian basin as Arctic analogues.
- Mean ocean floor age 64.2 Myr; spreading classes — Seton et al. (2020), G³ 21, e2020GC009214.
- Oceanic plateaus: 184, ~18.5 M km², ~5 % of the oceans — compilation summarised in "Oceanic plateau" (Wikipedia);
  Ontong Java 1.86 M km² — Neal, Coffin & Sager (2019), Oceanography 32(1).
- Wilson cycle, Iapetus (opened ~590 Ma, closed by Caledonian–Acadian orogenies ~440–400 Ma) — Dewey (2016) Geoscience Canada 43, 283; Wikipedia, "Iapetus Ocean".
- Old orogens keep relief for hundreds of Myr — Egholm et al. (2013) Nature 498, 475; Wolf et al. (2022) Nature 606, 516.
