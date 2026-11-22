# E01 — Global relief (bedrock elevation)

Script `scripts/relief.py`, parameters `data/relief.json`, references `research/R01_earth_references.md`. Input: the
tectonics layer (`tectonics.npz`, `boundaries.npz`, `volcanoes.npz`). Run: `python3 scripts/relief.py WORK_DIR`
(~2 min). Output `relief.npz` (6 × 1000 × 1000: `elevation_m`, `subgrid_relief_m`, `sediment_m`), not stored.

The first world built its relief from province levels and hand-chosen Earth analogues for each range. Here the
tectonics layer carries crustal thickness, so the mean elevation already comes from isostasy, and the relief layer
adds what isostasy does not give.

## Method

- **Continents:**
  - the isostatic surface (Airy on crustal thickness, swells, sea level);
  - a plain texture: a fractal field (wavelengths 31–2000 km) with a standard deviation by province, half the plain
    amplitudes of the first world (shields 125 m, platforms 60 m, basins 25 m, plateaus 150 m);
  - ridges and valleys: a ridged fractal (25–400 km) whose standard deviation is a share of the cell's elevation:
    young belts and arcs 0.4, old belts 0.25, plateaus 0.15 (smooth on top, like Tibet);
  - the ridges run along the belts. In fold-and-thrust belts (Appalachians, Zagros, Jura) ridges and valleys follow
    the strike and are several times longer than wide.
    - The across-strike direction is the gradient of crustal thickness, smoothed over ~50 km.
    - Each plane wave of the fractal is weighted at every cell by |cos|⁸ of the angle between its wave vector and
      that direction, so the crests run along strike.
    - The anisotropy is full where thickness changes by ≥ 20 m per km (belt flanks) and fades to isotropic on
      plateau tops and plains.
    - 20 % of the variance stays isotropic, which breaks the ridges into en-echelon segments.
- **Continental margins.** Where flooded continent lies within 400 km of ocean crust, it gets a shelf down to 130 m
  over the inner 40 % of the margin, then a smooth slope to the ocean floor, as in the first world. Wider flooded
  continent keeps its isostatic depth: these are inland seas.
- **Ocean floor:** the isostatic floor (GDH1 depth, thickened crust, swells), plus:
  - sediment (fit to GlobSed, R01), raising the floor by 0.57 of its thickness;
  - trenches at today's subduction zones: 3 km × exp(−(d/50 km)²), smoothed along strike;
  - residual depth (rms 261 m, 1000–4000 km) and abyssal hills (rms 75 m), as in the first world.
- **Volcanoes:** cones of radius r merge into the surface (taking the higher value). Each stands on the surface
  before any cone, so overlapping cones do not stack.
  - Hotspot centres on ocean crust: summit 4200 − 1300 √age m (primary) or 2800 − 1300 √age m (secondary), fitted to
    the Hawaiian chain. A volcano stays at least 3 km (1.5 km) above the floor, as the Emperor seamounts do.
    r = 60 km or 25 km.
  - Hotspot centres on continents: younger than 5 Myr only, +1500 m, r = 20 km.
  - Arc volcanoes: island arcs reach 0.5–2.5 km, r = 20 km; continental arcs stand 1–2.5 km above the arc,
    r = 15 km.
  - A cell mean is taken at no less than 4 km from the summit.
- **Relief inside a cell:** 40 m + a share of the elevation by province (young belts and arcs 0.5, old belts 0.35,
  plateaus 0.25, plains 0.05–0.1); on volcano flanks the slope × 10 km; 300 m on the sea floor.

## Runs

| Run | Change | Result |
|---|---|---|
| a | first version | volcano cones stacked on their neighbours (highest cell 25 km); a straight chain of islands 2500 km long (summits fell too slowly with age); trenches as rows of beads; a flood basalt as a flat block with straight scarps |
| b | cones stand on the surface before any cone; summits fitted to the Hawaiian chain (sea level at ~10–12 Ma) with a floor of 3 km above the sea floor; trenches smoothed along strike; flood basalts with irregular outlines and tapered rims (tectonics E05); continental arc volcanoes 1–2.5 km above the arc | ridges isotropic: no grain along the belts |
| c | ridges along strike (|cos|¹² weighting) | the grain follows the belts, but fine parallel stripes look combed |
| d | exponent 8 and a 20 % isotropic share | accepted: ridges along the belts, broken into segments (close-ups in `E01_relief_zoom.png`) |

## Result (run d)

| Check | This world | Earth |
|---|---|---|
| Land | 23.6 % of the surface | 29 % |
| Mean land elevation | 625 m | 840 m (with the ice sheets) |
| Land above 1000 / 2000 / 3000 m | 20 / 9.6 / 3.2 % of the land | ~27 / ~12 % (with the ice sheets) |
| Highest cell mean | 5829 m | 10 km means reach ~6 km in the Himalaya and Andes |
| Mean ocean depth | 3408 m | 3688 m |
| Ocean 0–500 / 500–2500 / 2500–4000 / 4000–5000 / 5000–6000 / >6000 m | 16.4 / 5.8 / 31.5 / 33.1 / 12.6 / 0.7 % | 10.2 / 9.2 / 27.5 / 33.2 / 19.2 / ~0.6 % |
| Shelf (sea floor above 200 m) | 29.1 M km² | 28.4 M km² |
| Deepest cell | −8386 m (trench) | −10 924 m (Challenger Deep) |
| Sediment: margins / far ocean / total | 3083 m / 583 m / 349 M km³ | 3044 m / 404 m / 337 M km³ |
| Volcanic islands | 90 hotspot, 331 arc; 0.13 M km² of land on ocean crust | — |

- **The ocean is ~280 m shallower than Earth's,** with more water over flooded continent (0–500 m: 16 % against
  10 %) and less deep floor. The ocean floor is young (mean 50 Myr against 64), and the sea stands 60 m higher.
  This is the Cretaceous-like state of a world 100 Myr after a breakup, not a missing feature.
- **The land is lower than Earth's.** The bedrock surface has no ice sheets, and the young orogens cover 8 % of the
  continental crust.
- **Shelf area and sediment** match Earth's within a few per cent, as calibrated in the first world.

## Known limits

- No river valleys, lakes or glacial erosion: the hydrology and ice layers will cut and fill the surface.
- Ridged noise dips below sea level in places inside low belts: small closed depressions, which the hydrology layer
  will turn into lakes or fill.
- A large primary swell on a continent makes a round dome ~1000 km across, as Earth's plume swells do (Hoggar,
  Tibesti, the East African plateau: ~1000 km across and ~1 km high).
- The continental arc volcanoes are spread evenly along the fronts (one per 70 km).
