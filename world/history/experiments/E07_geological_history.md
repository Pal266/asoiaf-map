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

Spreading half-rates: narrow sea 10 km/Myr, Mid-Sunset 55 (East Pacific Rise is faster, Mid-Atlantic
slower), South Sunset 15, Sothoryos–Ulthos 25, Antarctic ring 20.

## Runs

| Run | Change | Continental crust | Mean / max floor age, Ma | Mean basement depth, m | Verdict |
|---|---|---|---|---|---|
| a | first version, Parsons–Sclater depth | 34.4% | 84 / 180 | 5153 | no: each cell took the **youngest** capped age of any ridge on its plate, so distant floor got the narrow sea's 25 Ma; straight block outlines; unknown continents one uniform craton |
| b | nearest ridge per plate; fractal warp; craton cores and belts; shelves 200/60 km | 37.2% | 90 / 180 | 5259 | partly: shapes right; floor too old and deep because the Mid-Sunset ridge at 40 km/Myr needs 250 Myr to reach the far trenches (capped at 180 over a large area) |
| c | Mid-Sunset 55 km/Myr; GDH1 depth; passive margin 250 km; northern shield widened west | **38.7%** | **80 / 180** | **5023** | **yes** |

Earth reference: continental crust ~41% of the surface, oldest ocean floor ~180 Ma, typical ocean floor
~4.9 km below sea level. Run c is within a few percent; the remaining gap in depth is expected because
sediment and seamounts are not yet added. All 30 mapped canon features have the province E05 gives them.
Narrow sea floor is at most 25 Ma old.

## Known limits

- This is a scripted history, not a simulated one: no plate motion is integrated, so the ages of continental
  provinces are assigned by block, not derived.
- Age jumps across plate boundaries are sharp; on transforms that is real (fracture zones), elsewhere it is
  an artefact of the schematic E06 boundaries.
- GDH1 coefficients are the commonly quoted fitting equations; the paper itself was not opened in this session.
- Hot spots, oceanic plateaus and submerged microcontinents are not modelled.

## Conclusion

Run c is accepted as the history layer. It feeds the relief layer next: crust type and floor age give the
base elevation, province and event age give how high and rugged the land can be.

## Sources

- Parsons, B. & Sclater, J. G. (1977) JGR 82, 803 — age–depth relation.
- Stein, C. A. & Stein, S. (1992) Nature 359, 123 — GDH1 plate model.
- Oldest ocean floor ~180 Ma; spreading rates 2–15 cm/yr — physical geology textbooks (OpenText BC, Physical Geology 2e).
- Continental crust ~41% of Earth's surface — Wikipedia, "Continental crust".
- Wilson cycle, Iapetus (opened ~590 Ma, closed by Caledonian–Acadian orogenies ~440–400 Ma) — Dewey (2016) Geoscience Canada 43, 283; Wikipedia, "Iapetus Ocean".
- Old orogens keep relief for hundreds of Myr — Egholm et al. (2013) Nature 498, 475; Wolf et al. (2022) Nature 606, 516.
