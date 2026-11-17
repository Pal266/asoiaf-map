"""Land/sea mask from the official map "The Known World" (The Lands of Ice and Fire).

The map itself is not part of the repository (it is copyrighted); pass its path with --map.

Steps:
  1. Classify each pixel by colour: confident sea is blue (hue 190-255 deg, saturation >= 0.30, value >= 0.12).
     The pale shallow water drawn along the coasts (blue, saturation >= 0.10, value >= 0.62) is added where it
     connects to confident sea; the dark coastline stroke stops the growth. Everything else is land.
  2. Mask out what is not geography: the parchment frame and the title/legend box.
  2b. The Wall is drawn as a blue band of ice; a corridor along it is forced to land.
  3. Remove sea ice: the floes are the components lying wholly within a band along the top edge.
  4. Remove lettering: sea names are cream components with strokes thinner than ~13 px. Then fold lines: a
     morphological opening of the land and of the sea; land specks and sea holes smaller than a minimum area are
     absorbed by their surroundings.
  5. Measure the Wall between two control points given in map pixels and derive the scale.

Outputs (in --out): land_mask.png (white = land), preview.png, stats.json.
The same inputs always give the same outputs (no randomness).
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

Image.MAX_IMAGE_PIXELS = None

# Map pixels of the official 9050 x 6000 scan.
FRAME = 110                                  # parchment border width, px
LEGEND_BOX = (7090, 200, 8950, 1180)         # x0, y0, x1, y1 of the title and legend box
WALL_WEST = (1568, 1015)                     # Westwatch, where the Wall meets the Bay of Ice
WALL_EAST = (1933, 962)                      # Eastwatch, where the Wall meets the Bay of Seals
WALL_HALF_WIDTH = 18                         # px; the drawn band of ice is ~17 px wide
WALL_KM = 480.0                              # 100 leagues = 300 miles (canon, grade B)

SEA_HUE = (190.0, 255.0)
SEA_MIN_SAT = 0.30
SEA_MIN_VAL = 0.12
SHALLOW_MIN_SAT = 0.10                       # pale coastal water: light and weakly blue ...
SHALLOW_MIN_VAL = 0.62                       # ... and only where it connects to confident sea
OPEN_RADIUS = 4                              # px; removes strokes thinner than ~9 px (lettering, fold lines)
LABEL_ERODE = 6                              # px; lettering strokes (almost) vanish under this erosion
LABEL_CORE_FRAC = 0.12                       # a component keeping less than this share after erosion is lettering
MIN_ISLAND_PX = 150                          # land specks smaller than this become sea
MIN_LAKE_PX = 150                            # sea holes smaller than this become land
ICE_TOP_BAND = 900                           # px from the top in which white components are sea ice
ICE_MAX_PX = 200_000                         # larger white components are land (the Land of Always Winter)


def hsv(rgb):
    rgb = rgb.astype(np.float32) / 255.0
    mx = rgb.max(-1)
    mn = rgb.min(-1)
    d = mx - mn
    h = np.zeros_like(mx)
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    nz = d > 1e-6
    m = nz & (mx == r)
    h[m] = ((g - b)[m] / d[m]) % 6
    m = nz & (mx == g)
    h[m] = (b - r)[m] / d[m] + 2
    m = nz & (mx == b)
    h[m] = (r - g)[m] / d[m] + 4
    h *= 60.0
    s = np.where(mx > 0, d / np.maximum(mx, 1e-6), 0)
    return h, s, mx


def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def remove_small(mask, min_px):
    lab, n = ndimage.label(mask)
    if n == 0:
        return mask
    sizes = ndimage.sum(mask, lab, range(1, n + 1))
    keep = np.zeros(n + 1, bool)
    keep[1:] = sizes >= min_px
    return keep[lab]


def remove_lettering(land, h, s, v):
    lab, n = ndimage.label(land)
    ids = np.arange(1, n + 1)
    sizes = ndimage.sum(land, lab, ids)
    core = ndimage.binary_erosion(land, structure=disk(LABEL_ERODE))
    core_px = ndimage.sum(core, lab, ids)
    cream = land & (h >= 15) & (h <= 55) & (s < 0.5) & (v > 0.55)
    cream_frac = ndimage.sum(cream, lab, ids) / np.maximum(sizes, 1)
    label_ids = ids[(core_px / np.maximum(sizes, 1) < LABEL_CORE_FRAC) & (cream_frac > 0.5) & (sizes < 20_000)]
    return land & ~np.isin(lab, label_ids), int(len(label_ids))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True, help="path to the 9050 x 6000 scan of The Known World")
    ap.add_argument("--out", required=True, help="output directory")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rgb = np.asarray(Image.open(args.map).convert("RGB"))
    H, W = rgb.shape[:2]
    if (W, H) != (9050, 6000):
        raise SystemExit(f"expected a 9050 x 6000 scan, got {W} x {H}")
    h, s, v = hsv(rgb)

    # 1. colour classification: confident sea, then the pale shallow water along the coasts, grown from it
    blue = (h >= SEA_HUE[0]) & (h <= SEA_HUE[1])
    sea = blue & (s >= SEA_MIN_SAT) & (v >= SEA_MIN_VAL)
    shallow = blue & (s >= SHALLOW_MIN_SAT) & (v >= SHALLOW_MIN_VAL)
    sea = ndimage.binary_propagation(sea, mask=sea | shallow)

    # 2. not geography
    valid = np.zeros((H, W), bool)
    valid[FRAME:H - FRAME, FRAME:W - FRAME] = True
    x0, y0, x1, y1 = LEGEND_BOX
    valid[y0:y1, x0:x1] = False

    land = ~sea & valid

    # 2b. the Wall is drawn as a blue band of ice from coast to coast; it is land
    yy, xx = np.mgrid[min(WALL_WEST[1], WALL_EAST[1]) - 60:max(WALL_WEST[1], WALL_EAST[1]) + 60,
                      WALL_WEST[0] - 20:WALL_EAST[0] + 20]
    (ax, ay), (bx, by) = WALL_WEST, WALL_EAST
    t = np.clip(((xx - ax) * (bx - ax) + (yy - ay) * (by - ay)) / ((bx - ax) ** 2 + (by - ay) ** 2), 0, 1)
    dist = np.hypot(xx - (ax + t * (bx - ax)), yy - (ay + t * (by - ay)))
    wall = np.zeros((H, W), bool)
    wall[yy[dist <= WALL_HALF_WIDTH], xx[dist <= WALL_HALF_WIDTH]] = True
    land |= wall & valid

    # 3. sea ice: the floes are the components lying wholly in a band along the top edge (the large white
    #    Land of Always Winter reaches far below it)
    lab, n = ndimage.label(land)
    sizes = ndimage.sum(land, lab, range(1, n + 1))
    bottom = ndimage.maximum(np.indices((H, W))[0], lab, range(1, n + 1))
    ice_ids = np.where((bottom < ICE_TOP_BAND) & (sizes < ICE_MAX_PX))[0] + 1
    ice = np.isin(lab, ice_ids)
    land &= ~ice

    # 4. lettering: thin, uniformly cream components (sea names are written in cream on the sea)
    land, n_labels = remove_lettering(land, h, s, v)

    # 4b. fold lines and remaining specks
    land = ndimage.binary_opening(land, structure=disk(OPEN_RADIUS))
    sea2 = ndimage.binary_opening(~land & valid, structure=disk(OPEN_RADIUS))
    land = ~sea2 & valid
    land = remove_small(land, MIN_ISLAND_PX)
    holes = ~land & valid
    small_holes = holes & ~remove_small(holes, MIN_LAKE_PX)
    land |= small_holes
    land, n2 = remove_lettering(land, h, s, v)      # the opening thins some letters enough to show them
    n_labels += n2

    # 5. scale from the Wall
    wall_px = math.dist(WALL_WEST, WALL_EAST)
    km_per_px = WALL_KM / wall_px

    lab, n = ndimage.label(land)
    sizes = np.sort(ndimage.sum(land, lab, range(1, n + 1)))[::-1]
    stats = {
        "map_size_px": [W, H],
        "wall_px": round(wall_px, 1),
        "km_per_px_at_wall": round(km_per_px, 4),
        "map_extent_km_at_wall_scale": [round(W * km_per_px), round(H * km_per_px)],
        "land_fraction_of_valid_area": round(float(land.sum() / valid.sum()), 4),
        "land_components": int(n),
        "largest_components_km2_at_wall_scale": [round(float(a) * km_per_px ** 2) for a in sizes[:12]],
        "sea_ice_components_removed": int(len(ice_ids)),
        "label_components_removed": n_labels,
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    Image.fromarray((land * 255).astype(np.uint8)).save(out / "land_mask.png")

    prev = rgb.copy()
    prev[land] = (prev[land] * 0.5 + np.array([255, 255, 255]) * 0.5).astype(np.uint8)
    prev[~land & valid] = (prev[~land & valid] * 0.4).astype(np.uint8)
    edge = land ^ ndimage.binary_erosion(land, iterations=2)
    prev[edge] = [255, 0, 0]
    Image.fromarray(prev).resize((W // 4, H // 4), Image.LANCZOS).save(out / "preview.png")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
