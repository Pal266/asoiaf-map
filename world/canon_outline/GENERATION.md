# Layer "canon outline": how it is generated

The land/sea outline of the official map of the known world, placed on the planet and sampled onto the global grid.
The layer itself is not stored in the repository (Rule 4); these steps generate it, and the same input always gives
the same output.

## What the layer is

| File | Content |
|---|---|
| `canon_outline.npz` → `land_frac` | float32 `[6, 1000, 1000]`: share of the part of each cell covered by the map that is land, 0–1 |
| `canon_outline.npz` → `known_frac` | float32 `[6, 1000, 1000]`: share of each cell covered by the map, 0–1 (0 = outside the known world) |
| `islands.geojson` | the 286 islands smaller than 2 500 km² as polygons (longitude/latitude, degrees), with area, centroid and the grid cell of the centroid |
| `validation.json` | the checks below, with their results |
| `preview.png` | equirectangular preview, 3 600 × 1 800 px: sea blue, land beige, outside the map grey, 30° graticule |

The grid is the global cubed sphere of `world/grid/GRID.md` (6 faces × 1000 × 1000 cells, ~10 km, arrays indexed
`[face, j, i]`; face layout in `world/grid/scripts/cubesphere.py`). The map covers 16 % of the planet; elsewhere
`known_frac` is 0 and the outline will come from the geological history.

## Environment

Python 3.11 or newer (tested on 3.13) with the packages of `requirements.txt` at the repository root
(NumPy, SciPy, Pillow, scikit-image):

```
pip install -r requirements.txt
```

## Input

Map M1 (`research/SOURCES.md`): the official "The Known World", a 9050 × 6000 px JPEG scan of 28 766 695 bytes. It is
copyrighted and not part of the repository; keep it anywhere and pass its path. The scripts refuse any other size.

## Steps

All commands run from `world/canon_outline/scripts/`. `<work>` is any folder outside the repository.

### 1. Land/sea mask of the map — `land_mask.py`

```
python3 land_mask.py --map <path to M1> --out <work>/mask
```

About 40 s. Writes `land_mask.png` (9050 × 6000, white = land), `islands_px.json` (island polygons in map
pixels), `preview.png` and `stats.json`.

1. Each pixel is classified by colour. Confident sea: hue 190–255°, saturation ≥ 0.30, value ≥ 0.12. The pale
   shallow water along the coasts (blue, saturation ≥ 0.10, value ≥ 0.62) is added only where it connects to
   confident sea; the dark coastline stroke stops it. All other pixels are land.
2. The parchment frame (110 px) and the title and legend box are excluded (no data).
3. The Wall, drawn as a blue band of ice, is forced to land along the line between its two ends (half-width 18 px).
4. Sea ice: the land components lying wholly in the top 900 px (the floes of the Shivering Sea) become sea.
5. Lettering: sea names are cream components whose strokes nearly vanish under a 6 px erosion; they become sea.
6. Fold lines and specks: a morphological opening (radius 4 px) of land and of sea; land specks and sea holes below
   150 px are absorbed; the lettering test is repeated.
7. The scale: the Wall, from Westwatch (1568, 1015) to Eastwatch (1933, 962), 368.8 px for 480 km → 1.3014 km/px.
8. Islands. Land masses of at least 20 000 px keep the outline of steps 1–6. Every smaller land body is rebuilt
   from the classification before step 6:
   1. label ink (dark purple-maroon: value < 0.45, hue ≥ 255° or ≤ 15°), label halos (hue 190–260°, saturation
      < 0.22, value ≥ 0.6) and weakly saturated sea (hue 190–255°, saturation ≥ 0.25) are removed; ink is kept in
      the Smoking Sea box (3450–4400, 4300–5150 px), where the land itself is dark purple;
   2. a body is kept if, one pixel inside its edge, it holds at least 4 px of land-fill colours (khaki, green, tan,
      grey rock, snow, grey-blue rock shading; dark purple in the Smoking Sea box) and they are at least 25 % of its
      area;
   3. lettering is removed: bodies that almost vanish under a 6 px erosion and are more than 40 % pale (cream, pink,
      or cream blended into light sea), unless the darkest fifth of their edge is at least 0.07 darker than the sea
      around them (a coastline stroke);
   4. islands listed in `MANUAL_ISLANDS` are added from a box centred on each (now: Old Wyk);
   5. the bodies are closed (radius 2 px), their holes filled, and grown by 1 px to restore the coastline stroke.
9. Island polygons: every land body smaller than 2 500 km² is traced (marching squares at the 0.5 level, simplified
   to 0.5 px).

Why each threshold has its value: experiments `experiments/E01_land_mask.md` (steps 1–7) and
`experiments/E03_small_islands.md` (steps 8–9).

### 2. The outline on the global grid — `canon_outline.py`

```
python3 canon_outline.py --mask <work>/mask/land_mask.png --islands <work>/mask/islands_px.json \
                         --out <work>/canon_outline
```

About 45 s.

1. Placement (decision 3 of `world/canon/CONSTRAINTS.md`, experiment `experiments/E02_projection.md`): the map is a
   sinusoidal projection centred on the Wall. The Wall's midpoint (1750.5, 988.5 px) lies at 63.75° N on the
   meridian 0° (the longitude is a free choice). For a point at latitude φ and longitude λ:
   `x = 1750.5 + R·λ·cos φ / 1.3014`, `y = 988.5 + R·(63.75° − φ) / 1.3014`, with R = 6 371 km and angles in
   radians.
2. North of ~62° N the map is wider than the planet's circumference, so its north-east corner lies beyond the
   antimeridian: each point is looked up at λ, then λ + 360°, then λ − 360°, and the first position on the map
   with data is used.
3. Each cell is sampled at 8 × 8 points spread evenly in its equiangular coordinates (~1.2 km apart, about one map
   pixel). `known_frac` is the share of points that fall on the map with data; `land_frac` is the share of those
   that are land.
4. The island polygons are converted from map pixels to longitude/latitude with the same placement.

## Checks (in `validation.json`)

| Check | Expected | Result |
|---|---|---|
| Grid area = sphere area | equal | 510 064 472 km² both |
| Area covered by the map, grid vs map (the projection is equal-area) | equal | 83.35 vs 83.37 M km², 0.02 % (the map overlaps itself north of 73.3° N, over open sea) |
| Land area, grid vs map | equal | 32.1035 vs 32.1034 M km², < 0.001 % |
| Six largest land masses, grid vs map | equal | Essos part 21.37 M km², Westeros 7.04 M km², then 1.27, 0.50, 0.42, … M km²; all within 0.01 % |
| Islands below 2 500 km² | stored as polygons | 286 islands, 161 000 km²; 23 smaller than one cell |
| The Wall's length on the globe | 480 km (canon) | 479.6 km |
| Seams between cube faces | no jump | mean land-share difference 0.0245 across face edges vs 0.0249 one row inside |
| Determinism | identical output | two runs of each step give byte-identical outputs |

Latitudes and longitudes of places read on M1: Castle Black 63.8° N 1.3° E, Winterfell 56.5° N 4.3° W,
King's Landing 38.7° N 1.2° E, Oldtown 30.0° N 11.1° W, Sunspear 28.7° N 5.0° E, Volantis 28.7° N 22.9° E,
Qarth 22.5° N 57.3° E, Asshai 9.0° N 83.0° E.

## Known limits

- Islets drawn as dark specks without a land-coloured fill are lost (many in the Thousand Islands, mostly under
  ~40 km²); the smallest island kept is 44 km² (E03).
- A map label lying over the coast of a large land mass stays part of that coast (E01, E03).
- Lakes are not separated from the sea: dark blue lakes count as sea, small ones as land. Lakes need their own layer.
- Some noise remains in the Land of Always Winter (E01).
