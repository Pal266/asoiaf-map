# E04. Canon relief from the official map

## Goal

Record where canon puts mountains, hills, volcanoes, deserts, basins, inland lakes and wetlands, so that the geological
history can be built to explain them (stage 1 of the geological history).

## Input data

- Map M1 (`world/canon_outline/research/SOURCES.md`) and the land mask of the canon outline layer (E03).
- Canon facts on A Wiki of Ice and Fire (each feature's `source` in `data/canon_features.json`).

## Method and iterations

| Run | Method | Result |
|---|---|---|
| a | Fully automatic: density of relief strokes (pixels at least 0.12 darker or 0.10 brighter than the local background), forest and desert by colour | Mountains found (Bone Mountains, Frostfangs, Mountains of the Moon), but coasts, rivers, roads, labels, town icons, the fold lines and forests also produced strokes |
| b | Strokes excluded on lines (rivers, roads, labels, icons, dilated 3 px), a 10 px coastal band, the fold lines (rows 1985 and 3986, columns 1650 and 7240, found as bright lines across the sea) and forests | Cleaner, but town icons and halos still read as hills, and relief drawn in brown on orange land (the Red Mountains) or in light strokes (the Painted Mountains) was weak |
| c | **Envelopes drawn by eye** on gridded views of M1 for every named or clearly drawn feature (34 features), with canon facts and sources; the stroke density of run b kept only inside mountain, hill and volcanic envelopes, as a measure of how strongly relief is drawn | Final run |

The automatic density alone cannot tell relief from the map's other drawings, while the features are few and large;
envelopes by eye are reliable at the scale of the geological history, and the density inside them keeps the drawn
detail for the later relief layer.

## Results (run c)

34 features: 10 mountain ranges or massifs (one a point: the Mother of Mountains), 12 hill regions, 3 volcanic (the
Fourteen Flames and the Valyrian peninsula, Dragonmont, Marahai), 4 deserts, 1 canyon basin (the Dry Deep), 2 inland
lakes (the Shrinking Sea, the Bleeding Sea), 1 wetland (the Neck) and the Stepstones as the remnant of the Arm of
Dorne.

| Type | Share of the land on the map |
|---|---|
| Mountain envelopes | 16.2 % |
| Hill envelopes | 7.0 % |
| Desert envelopes | 7.6 % |
| Volcanic | 1.8 % |
| Basin, wetland | 1.5 % |
| Strongly drawn relief (density > 0.3) | 1.5 % |

Latitudes on the globe (decision 3 of `world/canon/CONSTRAINTS.md`): the Frostfangs 64–73° N, the Mountains of the
Moon 43–50° N, the Red Mountains 29–36° N, the Bone Mountains 24–49° N, the deserts of Dorne, the red waste and the
Great Sand Sea 23–40° N, the Shadow Lands 10–20° N.

Grid areas match the map areas within 3 % for every feature (most within 0.5 %); two runs give identical outputs.

## Conclusions

- **Goal reached.** The canon relief is recorded as envelopes and points in longitude/latitude and on the global grid.
- Mountain envelopes cover 16 % of the known land, hills 7 %; Earth's mountain share is ~24 % by the UNEP-WCMC
  definition (Kapos et al. 2000), so the known world is somewhat less mountainous than Earth, within a plausible range.
- The relief density is a weak indicator in brown-on-orange areas (Red Mountains, 0.04) and must not be read as
  height; heights come from the geological history.
- Not recorded here: forests (a vegetation target, later), rivers and lakes other than the two inland seas (a
  hydrology target, later), and relief too faint to see on M1.
