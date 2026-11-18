"""Layer "canon relief": the relief features drawn on the official map, as targets for the geological history.

The process is described in world/canon_relief/GENERATION.md.

Input:
  --map       map M1 (9050 x 6000 px), not part of the repository
  --mask      land_mask.png from world/canon_outline/scripts/land_mask.py
  data/canon_features.json (in the repository): named envelopes and points drawn by eye on M1
Output (in --out):
  canon_relief.npz     feature_id int16 [6, 1000, 1000]: id of the canon feature at the cell centre (0 = none)
                       relief    float32 [6, 1000, 1000]: density of drawn relief strokes inside mountain and hill
                                 envelopes, 0..1 (mean of 4 x 4 samples per cell; 0 elsewhere)
  canon_features.geojson   envelopes (polygons) and points in longitude/latitude, with name, type, area, latitudes
  stats.json           per-feature measures and checks
  preview_map.png      envelopes and relief drawn over M1 (2262 x 1500)
Deterministic: the same inputs give the same files.
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from skimage.draw import polygon as fill_polygon

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "canon_outline" / "scripts"))
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
import cubesphere as cs                                       # noqa: E402
import land_mask as lm                                        # noqa: E402
from canon_outline import KM_PER_PX, globe, map_px            # noqa: E402

Image.MAX_IMAGE_PIXELS = None
FEATURES = HERE.parent / "data" / "canon_features.json"
CLIP_TO_LAND = {"mountains", "hills", "desert", "volcanic", "wetland", "basin"}
RELIEF_TYPES = {"mountains", "hills", "volcanic"}
SUB = 4

# Drawn relief strokes (experiment E04): shadows and highlights against the local background, excluding lines
# (rivers, roads, labels, icons), the coastal band, the map's fold lines and forests.
BG_SIGMA = 15
SHADOW_DARKER = 0.12
HIGHLIGHT_BRIGHTER = 0.10
DENSITY_SIGMA = 15
FOLD_ROWS = (1985, 3986)
FOLD_COLS = (1650, 7240)


def relief_density(rgb, land):
    h, s, v = lm.hsv(rgb)
    bg = ndimage.gaussian_filter(v, BG_SIGMA)
    lines = (((h >= 185) & (h <= 265) & (s >= 0.12) & (v < 0.7))
             | (((h >= 290) | (h <= 30)) & (s >= 0.2) & (v < 0.65)))
    lines = ndimage.binary_dilation(lines, iterations=3)
    forest = ndimage.gaussian_filter(((h >= 60) & (h <= 170) & (s >= 0.12) & (v < 0.5)).astype(np.float32), 10) > 0.3
    coast = ndimage.binary_dilation(~land, iterations=10)
    fold = np.zeros_like(land)
    for y in FOLD_ROWS:
        fold[y - 20:y + 21] = True
    for x in FOLD_COLS:
        fold[:, x - 12:x + 13] = True
    usable = land & ~(lines | coast | fold | forest)
    strokes = usable & ((v < bg - SHADOW_DARKER) | ((v > bg + HIGHLIGHT_BRIGHTER) & (s < 0.25)))
    num = ndimage.gaussian_filter(strokes.astype(np.float32), DENSITY_SIGMA)
    den = ndimage.gaussian_filter(usable.astype(np.float32), DENSITY_SIGMA)
    return np.where(den > 0.5, num / np.maximum(den, 1e-6), 0).astype(np.float32)


def lookup(lat, lon, valid_shape):
    """Map pixel of globe points, trying lon, lon + 360, lon - 360 (see canon_outline.map_px)."""
    H, W = valid_shape
    xi = np.full(lat.shape, -1, np.int64)
    yi = np.full(lat.shape, -1, np.int64)
    found = np.zeros(lat.shape, bool)
    for turn in (0, 1, -1):
        x, y = map_px(lat, lon, turn)
        xt, yt = np.floor(x).astype(np.int64), np.floor(y).astype(np.int64)
        ok = ~found & (xt >= 0) & (xt < W) & (yt >= 0) & (yt < H)
        xi[ok], yi[ok] = xt[ok], yt[ok]
        found |= ok
    return xi, yi, found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rgb = lm.load_map(args.map)
    land = np.asarray(Image.open(args.mask)) > 0
    H, W = land.shape
    feats = json.loads(FEATURES.read_text())["features"]

    # ---- envelopes on the map ----
    fid = np.zeros((H, W), np.int16)
    for f in feats:
        if "polygon_px" not in f:
            continue
        p = np.array(f["polygon_px"], float)
        rr, cc = fill_polygon(p[:, 1], p[:, 0], shape=(H, W))
        m = np.zeros((H, W), bool)
        m[rr, cc] = True
        if f["type"] in CLIP_TO_LAND:
            m &= land
        fid[m & (fid == 0)] = f["id"]
    dens = relief_density(rgb, land)
    relief_types = np.array([0] + [f["id"] for f in feats if f["type"] in RELIEF_TYPES])
    in_relief = np.isin(fid, relief_types[1:])
    relief = np.where(in_relief, np.clip(dens, 0, 1), 0).astype(np.float32)

    # ---- per-feature measures ----
    px_km2 = KM_PER_PX ** 2
    stats_f = []
    for f in feats:
        rec = {"id": f["id"], "name": f["name"], "type": f["type"], "region": f["region"]}
        if "polygon_px" in f:
            m = fid == f["id"]
            ys, xs = np.nonzero(m)
            lat, _ = globe(xs[::50] + 0.5, ys[::50] + 0.5)
            rec["area_km2"] = round(float(m.sum()) * px_km2)
            rec["lat_range"] = [round(math.degrees(lat.min()), 1), round(math.degrees(lat.max()), 1)]
            if f["type"] in RELIEF_TYPES:
                rec["relief_mean"] = round(float(dens[m].mean()), 3)
                rec["relief_share_above_0.3"] = round(float((dens[m] > 0.3).mean()), 3)
        else:
            la, lo = globe(f["point_px"][0] + 0.5, f["point_px"][1] + 0.5)
            rec["lat_lon"] = [round(math.degrees(la), 2), round(math.degrees(lo), 2)]
        stats_f.append(rec)

    # ---- grid ----
    n = cs.N
    e = cs.edges(n)
    d = e[1] - e[0]
    grid_fid = np.zeros((6, n, n), np.int16)
    grid_relief = np.zeros((6, n, n), np.float32)
    for face in range(6):
        t, u = np.meshgrid(np.tan(e[:-1] + 0.5 * d), np.tan(e[:-1] + 0.5 * d))
        lat, lon = cs.xyz_to_latlon(*cs.face_xyz(face, t, u))
        xi, yi, ok = lookup(lat, lon, (H, W))
        g = np.zeros((n, n), np.int16)
        g[ok] = fid[yi[ok], xi[ok]]
        grid_fid[face] = g
        acc = np.zeros((n, n))
        for oj in (np.arange(SUB) + 0.5) / SUB:
            for oi in (np.arange(SUB) + 0.5) / SUB:
                t, u = np.meshgrid(np.tan(e[:-1] + oi * d), np.tan(e[:-1] + oj * d))
                lat, lon = cs.xyz_to_latlon(*cs.face_xyz(face, t, u))
                xi, yi, ok = lookup(lat, lon, (H, W))
                r = np.zeros((n, n))
                r[ok] = relief[yi[ok], xi[ok]]
                acc += r
        grid_relief[face] = acc / SUB ** 2
    np.savez_compressed(out / "canon_relief.npz", feature_id=grid_fid, relief=grid_relief)

    area = cs.cell_area_km2(n)
    A = np.broadcast_to(area, (6, n, n))
    grid_area = {}
    for f in feats:
        if "polygon_px" in f:
            grid_area[f["id"]] = float(A[grid_fid == f["id"]].sum())
    for rec in stats_f:
        if rec["id"] in grid_area:
            rec["area_km2_grid"] = round(grid_area[rec["id"]])

    land_km2 = float(land.sum()) * px_km2
    by_type = {}
    for f in feats:
        if "polygon_px" in f:
            by_type[f["type"]] = by_type.get(f["type"], 0) + float((fid == f["id"]).sum()) * px_km2
    stats = {
        "features": stats_f,
        "envelope_area_by_type_km2": {k: round(v) for k, v in sorted(by_type.items())},
        "land_km2_on_map": round(land_km2),
        "share_of_land": {k: round(v / land_km2, 4) for k, v in sorted(by_type.items()) if k in CLIP_TO_LAND},
        "drawn_relief_above_0.3_share_of_land": round(float(((dens > 0.3) & in_relief).sum()) * px_km2 / land_km2, 4),
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    # ---- lon/lat geometry ----
    gj = []
    for f, rec in zip(feats, stats_f):
        props = {k: v for k, v in rec.items()}
        props["canon"] = f["canon"]
        props["source"] = f["source"]
        if "polygon_px" in f:
            p = np.array(f["polygon_px"] + [f["polygon_px"][0]], float)
            la, lo = globe(p[:, 0], p[:, 1])
            lo = (lo + np.pi) % (2 * np.pi) - np.pi
            geom = {"type": "Polygon",
                    "coordinates": [[[round(math.degrees(b), 4), round(math.degrees(a), 4)] for a, b in zip(la, lo)]]}
        else:
            la, lo = globe(*f["point_px"])
            geom = {"type": "Point", "coordinates": [round(math.degrees(lo), 4), round(math.degrees(la), 4)]}
        gj.append({"type": "Feature", "properties": props, "geometry": geom})
    (out / "canon_features.geojson").write_text(
        json.dumps({"type": "FeatureCollection", "features": gj}, ensure_ascii=False))

    # ---- preview over the map ----
    colours = {"mountains": (255, 40, 40), "hills": (255, 200, 0), "volcanic": (255, 0, 255), "desert": (255, 140, 0),
               "basin": (160, 80, 0), "lake": (0, 200, 255), "wetland": (0, 255, 120), "sunken_bridge": (180, 180, 255)}
    img = (rgb * 0.6).astype(np.uint8)
    rmask = relief > 0.3
    img[rmask] = (img[rmask] * 0.4 + np.array([255, 60, 60]) * 0.6).astype(np.uint8)
    prev = Image.fromarray(img).resize((2262, 1500), Image.LANCZOS)
    dr = ImageDraw.Draw(prev)
    sc = 2262 / W
    for f in feats:
        c = colours[f["type"]]
        if "polygon_px" in f:
            pts = [(x * sc, y * sc) for x, y in f["polygon_px"]]
            dr.line(pts + [pts[0]], fill=c, width=2)
            cx = sum(p[0] for p in pts) / len(pts)
            cy = sum(p[1] for p in pts) / len(pts)
        else:
            cx, cy = f["point_px"][0] * sc, f["point_px"][1] * sc
            dr.ellipse((cx - 5, cy - 5, cx + 5, cy + 5), outline=c, width=2)
        dr.text((cx + 6, cy - 6), str(f["id"]), fill=c)
    prev.save(out / "preview_map.png")
    print(json.dumps({k: v for k, v in stats.items() if k != "features"}, indent=1))


if __name__ == "__main__":
    main()
