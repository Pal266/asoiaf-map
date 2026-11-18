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
  6. Restore islands (experiment E03). Steps 1-4 keep only land masses of at least BIG_LAND_PX; every smaller land
     body is rebuilt from the colour classification before step 4b (the opening that erased small islands):
     a. label ink (dark purple-maroon, the island and place names), label halos (light, weakly blue) and weakly
        saturated sea (blue, saturation >= 0.25) are removed; ink is kept in the Smoking Sea, where the land itself is
        dark purple;
     b. each remaining body is kept if its interior holds land-fill colours (khaki, green, tan, grey rock, snow;
        purple only in the Smoking Sea): at least ISLAND_MIN_FILL_PX pixels and ISLAND_MIN_FILL_SHARE of its area;
     c. lettering (thin and mostly pale: cream, pink, or cream blended into light sea) is removed, unless a darker
        outline surrounds the body (the darkest fifth of its edge at least ISLAND_OUTLINE_DARKER darker than the sea around it);
     d. the kept bodies are closed (radius 2 px), their holes filled, and grown by 1 px to restore the outline
        stroke that step a may have thinned.
     e. named islands the classifier misses are added from boxes listed in MANUAL_ISLANDS (Old Wyk).
  7. Island polygons: every land body smaller than ISLAND_VECTOR_MAX_KM2 is traced as a polygon (map px).

Outputs (in --out): land_mask.png (white = land), islands_px.json (island polygons in map px), preview.png,
stats.json.
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

# Step 6, islands (E03)
BIG_LAND_PX = 20_000                         # land masses at least this large come from steps 1-4b
ISLAND_MIN_FILL_PX = 4                       # land-fill pixels needed inside a small land body ...
ISLAND_MIN_FILL_SHARE = 0.25                 # ... and their share of its area
ISLAND_MAX_PX = 40_000
SEAISH_MIN_SAT = 0.25                        # blue at least this saturated is sea, not grey-blue rock shading
ISLAND_OUTLINE_DARKER = 0.07                 # a thin pale body with an edge this much darker than the sea is an island
SMOKING_SEA_BOX = (3450, 4300, 4400, 5150)   # x0, y0, x1, y1: dark purple Valyrian land, ink test not applied
ISLAND_VECTOR_MAX_KM2 = 2500.0
# Named islands that the classifier misses (checked by eye on M1, experiment E03). Each box is centred on the island;
# the island is the body at the centre that is neither confident sea nor label ink (after a 1 px opening).
MANUAL_ISLANDS = {
    "Old Wyk": (930, 2445, 990, 2505),
}               # bodies smaller than this are also stored as polygons


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


def land_fill(h, s, v):
    """Colours of land on M1: khaki/green/tan/orange, grey rock, snow, grey-blue rock shading.
    (Purple Valyrian land is added only inside SMOKING_SEA_BOX, see restore_islands: elsewhere mauve is lettering.)"""
    return (((h >= 20) & (h <= 130) & (s >= 0.2) & (v >= 0.35))
            | ((h >= 15) & (h <= 90) & (s < 0.25) & (v >= 0.45))
            | ((s < 0.08) & (v >= 0.8))
            | ((h >= 190) & (h <= 260) & (s < 0.2) & (v >= 0.4) & (v < 0.6)))


def restore_islands(land, raw_land, valid, h, s, v):
    """Step 6: rebuild every land body smaller than BIG_LAND_PX. Returns the new mask and counters."""
    H, W = land.shape
    lab, _ = ndimage.label(land)
    sz = np.bincount(lab.ravel())
    sz[0] = 0
    big = (sz >= BIG_LAND_PX)[lab]

    ink = (v < 0.45) & ((h >= 255) | (h <= 15))
    halo = (h >= 190) & (h <= 260) & (s < 0.22) & (v >= 0.6)
    x0, y0, x1, y1 = SMOKING_SEA_BOX
    valyria = np.zeros((H, W), bool)
    valyria[y0:y1, x0:x1] = True
    fill = land_fill(h, s, v) | (valyria & (h >= 250) & (v >= 0.2))
    ink &= ~valyria

    seaish = (h >= 190) & (h <= 255) & (s >= SEAISH_MIN_SAT)   # weakly saturated sea that step 1 left as land
    cand = raw_land & ~ink & ~halo & ~seaish & ~ndimage.binary_dilation(big, iterations=2)
    cl, cn = ndimage.label(cand)
    ids = np.arange(1, cn + 1)
    csz = ndimage.sum(cand, cl, ids)
    inner_fill = ndimage.binary_erosion(cand, iterations=1) & fill
    nf = ndimage.sum(inner_fill, cl, ids)
    keep = ids[(nf >= ISLAND_MIN_FILL_PX) & (nf >= ISLAND_MIN_FILL_SHARE * csz) & (csz < ISLAND_MAX_PX)]
    isl = np.isin(cl, keep)

    # lettering: thin and mostly pale (cream, pink, or cream blended into light sea); khaki islands are more saturated
    il, iln = ndimage.label(isl)
    iids = np.arange(1, iln + 1)
    isz = ndimage.sum(isl, il, iids)
    core = ndimage.binary_erosion(isl, structure=disk(LABEL_ERODE))
    core_frac = ndimage.sum(core, il, iids) / np.maximum(isz, 1)
    pale = (((s < 0.32) & (v > 0.6) & ~((h >= 190) & (h <= 255) & (s >= 0.2)))
            | ((s < 0.5) & (v > 0.55) & ((h <= 20) | (h >= 330))))      # pink lettering
    pale_frac = ndimage.sum(pale & isl, il, iids) / np.maximum(isz, 1)
    dark_frac = ndimage.sum((v < 0.45) & isl, il, iids) / np.maximum(isz, 1)
    suspects = iids[(core_frac < LABEL_CORE_FRAC) & (pale_frac > 0.4)]
    # a thin pale body is still an island if a darker outline surrounds it (the coastline stroke)
    objs = ndimage.find_objects(il)
    letters = []
    for c in suspects:
        sl = objs[c - 1]
        ya, yb = max(0, sl[0].start - 10), min(H, sl[0].stop + 10)
        xa, xb = max(0, sl[1].start - 10), min(W, sl[1].stop + 10)
        m = il[ya:yb, xa:xb] == c
        ring = ndimage.binary_dilation(m, iterations=2) & ~ndimage.binary_erosion(m, iterations=1)
        outer = ndimage.binary_dilation(m, structure=disk(8)) & ~ndimage.binary_dilation(m, structure=disk(4))
        vv = v[ya:yb, xa:xb]
        if not (outer.any() and np.percentile(vv[ring], 20) < np.median(vv[outer]) - ISLAND_OUTLINE_DARKER):
            letters.append(c)
    isl &= ~np.isin(il, letters)

    for name, (bx0, by0, bx1, by1) in MANUAL_ISLANDS.items():
        hb, sb, vb = h[by0:by1, bx0:bx1], s[by0:by1, bx0:bx1], v[by0:by1, bx0:bx1]
        body = ~((hb >= SEA_HUE[0]) & (hb <= SEA_HUE[1]) & (sb >= SEA_MIN_SAT)) & ~ink[by0:by1, bx0:bx1]
        bl, bn = ndimage.label(ndimage.binary_opening(body))
        centre = bl[(by1 - by0) // 2, (bx1 - bx0) // 2]
        if centre:
            isl[by0:by1, bx0:bx1] |= bl == centre

    isl = ndimage.binary_closing(isl, structure=disk(2))
    isl = ndimage.binary_fill_holes(isl)
    isl = ndimage.binary_dilation(isl, iterations=1) & valid
    return big | isl, {"island_bodies_kept": int(len(keep)), "island_lettering_removed": int(len(letters)),
                       "manual_islands": list(MANUAL_ISLANDS)}


def island_polygons(land, km_per_px):
    """Step 7: polygons (map px, x/y) of land bodies smaller than ISLAND_VECTOR_MAX_KM2."""
    from skimage import measure
    lab, n = ndimage.label(land)
    objs = ndimage.find_objects(lab)
    sizes = np.bincount(lab.ravel())
    out = []
    max_px = ISLAND_VECTOR_MAX_KM2 / km_per_px ** 2
    for i in range(1, n + 1):
        if sizes[i] >= max_px:
            continue
        sl = objs[i - 1]
        m = np.pad(lab[sl] == i, 1)
        contours = measure.find_contours(m.astype(float), 0.5)
        ring = max(contours, key=len)
        ring = measure.approximate_polygon(ring, tolerance=0.5)
        ys = ring[:, 0] - 1 + sl[0].start
        xs = ring[:, 1] - 1 + sl[1].start
        cy, cx = ndimage.center_of_mass(lab[sl] == i)
        out.append({"area_px": int(sizes[i]), "area_km2": round(float(sizes[i]) * km_per_px ** 2, 2),
                    "centroid_px": [round(cx + sl[1].start, 1), round(cy + sl[0].start, 1)],
                    "polygon_px": [[round(float(x), 1), round(float(y), 1)] for x, y in zip(xs, ys)]})
    out.sort(key=lambda d: (d["centroid_px"][1], d["centroid_px"][0]))
    for k, d in enumerate(out):
        d["id"] = k + 1
    return out


def load_map(path):
    rgb = np.asarray(Image.open(path).convert("RGB"))
    H, W = rgb.shape[:2]
    if (W, H) != (9050, 6000):
        raise SystemExit(f"expected a 9050 x 6000 scan, got {W} x {H}")
    return rgb


def build(rgb):
    """Steps 1-4: returns the final land mask, the mask before the opening (raw_land), validity and counters."""
    H, W = rgb.shape[:2]
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
    raw_land = land.copy()

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

    land_e01 = land
    land, isl_stats = restore_islands(land_e01, raw_land, valid, h, s, v)
    return {"land": land, "land_e01": land_e01, "raw_land": raw_land, "valid": valid, "ice_ids": ice_ids,
            "n_labels": n_labels, "islands": isl_stats, "hsv": (h, s, v)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True, help="path to the 9050 x 6000 scan of The Known World")
    ap.add_argument("--out", required=True, help="output directory")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rgb = load_map(args.map)
    H, W = rgb.shape[:2]
    r = build(rgb)
    land, valid, ice_ids, n_labels = r["land"], r["valid"], r["ice_ids"], r["n_labels"]

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
        **r["islands"],
    }
    islands = island_polygons(land, km_per_px)
    stats["island_polygons"] = len(islands)
    stats["island_polygons_area_km2"] = round(sum(d["area_km2"] for d in islands), 1)
    (out / "islands_px.json").write_text(json.dumps(islands))
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
