# Layer "canon outline": how it is generated

The land/sea outline of the official map of the known world, placed on the planet and sampled onto the global grid.
The layer itself is not stored in the repository (Rule 4); these steps generate it, and the same input always gives
the same output.

## What the layer is

| File | Content |
|---|---|
| `canon_outline.npz` → `land_frac` | float32 `[6, 1000, 1000]`: share of the part of each cell covered by the map that is land, 0–1 |
| `canon_outline.npz` → `known_frac` | float32 `[6, 1000, 1000]`: share of each cell covered by the map, 0–1 (0 = outside the known world) |
| `validation.json` | the checks below, with their results |
| `preview.png` | equirectangular preview, 3 600 × 1 800 px: sea blue, land beige, outside the map grey, 30° graticule |

The grid is the global cubed sphere of `world/grid/GRID.md` (6 faces × 1000 × 1000 cells, ~10 km, arrays indexed
`[face, j, i]`; face layout in `world/grid/scripts/cubesphere.py`). The map covers 16 % of the planet; elsewhere
`known_frac` is 0 and the outline will come from the geological history.

## Environment

Python 3.11 or newer (tested on 3.13) with the packages of `requirements.txt` at the repository root:

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

About 20 s. Writes `land_mask.png` (9050 × 6000, white = land), `preview.png` and `stats.json`.

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

Why each threshold has its value: experiment `experiments/E01_land_mask.md`.

### 2. The outline on the global grid — `canon_outline.py`

```
python3 canon_outline.py --mask <work>/mask/land_mask.png --out <work>/canon_outline
```

About 10 s.

1. Placement (decision 3 of `world/canon/CONSTRAINTS.md`, experiment `experiments/E02_projection.md`): the map is a
   sinusoidal projection centred on the Wall. The Wall's midpoint (1750.5, 988.5 px) lies at 63.75° N on the
   meridian 0° (the longitude is a free choice). For a point at latitude φ and longitude λ:
   `x = 1750.5 + R·λ·cos φ / 1.3014`, `y = 988.5 + R·(63.75° − φ) / 1.3014`, with R = 6 371 km and angles in
   radians.
2. North of ~62° N the map is wider than the planet's circumference, so its north-east corner lies beyond the
   antimeridian: each point is looked up at λ, then λ + 360°, then λ − 360°, and the first position on the map
   with data is used.
3. Each cell is sampled at 3 × 3 points spread evenly in its equiangular coordinates. `known_frac` is the share of
   points that fall on the map with data; `land_frac` is the share of those that are land.

## Checks (in `validation.json`)

| Check | Expected | Result |
|---|---|---|
| Grid area = sphere area | equal | 510 064 472 km² both |
| Area covered by the map, grid vs map (the projection is equal-area) | equal | 83.35 vs 83.37 M km², 0.02 % (the map overlaps itself north of 73.3° N, over open sea) |
| Land area, grid vs map | equal | 32.180 vs 32.179 M km², 0.003 % |
| Six largest land masses, grid vs map | equal | Essos part 21.37 M km², Westeros 7.04 M km², then 1.27, 0.50, 0.42, 0.26 M km²; all within 0.04 % |
| The Wall's length on the globe | 480 km (canon) | 479.6 km |
| Seams between cube faces | no jump | mean land-share difference 0.0233 across face edges vs 0.0234 one row inside |
| Determinism | identical output | two runs give byte-identical `canon_outline.npz` |

Latitudes and longitudes of places read on M1: Castle Black 63.8° N 1.3° E, Winterfell 56.5° N 4.3° W,
King's Landing 38.7° N 1.2° E, Oldtown 30.0° N 11.1° W, Sunspear 28.7° N 5.0° E, Volantis 28.7° N 22.9° E,
Qarth 22.5° N 57.3° E, Asshai 9.0° N 83.0° E.

## Known limits

- Islands smaller than ~250 km² on the map are lost (E01).
- Lakes are not separated from the sea: dark blue lakes count as sea, small ones as land. Lakes need their own layer.
- Some noise remains in the Thousand Islands and the Land of Always Winter (E01).
