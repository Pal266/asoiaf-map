"""The main rivers and inland lakes of the official map on the cubed sphere (experiment E12).

Each river of data/canon_rivers.json is a polyline traced by eye on M1. It is snapped to the river line drawn on the
map: every vertex moves to the darkest pixel within 12 px, then each segment follows the least-cost path through the
map's dark pixels (cost ~ brightness^3) inside a corridor of corridor_px around the traced line. The path pixels are
placed on the globe (world/canon_outline: sinusoidal projection centred on the Wall) and on the 10 km cells.

Input:
  --map   M1 (world/canon_outline/research/SOURCES.md)
Output (in --out):
  canon_rivers.npz  river_id int16 [6, 1000, 1000]: 1 + index of the river in canon_rivers.json, 0 elsewhere;
                    r<k>: flat cell indices of river k in order from source to mouth;
                    lake_id int16: 1 + index of the canon lake, 0 elsewhere
  canon_rivers.json the snapped paths in latitude/longitude (for checks), with lengths
  check.png         the snapped paths over M1 (crops of Westeros and western Essos)
Deterministic.
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from skimage.graph import route_through_array

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
sys.path.insert(0, str(WORLD / "canon_outline" / "scripts"))
import cubesphere as cs                    # noqa: E402
from canon_outline import KM_PER_PX, globe  # noqa: E402
from land_mask import hsv                   # noqa: E402

Image.MAX_IMAGE_PIXELS = None
DATA = HERE.parent / "data" / "canon_rivers.json"
SNAP_PX = 12


def snap_river(v, pts, corridor):
    """Least-cost path through dark pixels along a traced polyline. v: brightness (0-1) of the whole map."""
    pts = np.array(pts, float)
    x0, y0 = (pts.min(0) - corridor - SNAP_PX - 2).astype(int)
    x1, y1 = (pts.max(0) + corridor + SNAP_PX + 2).astype(int)
    x0, y0 = max(x0, 0), max(y0, 0)
    sub = v[y0:y1, x0:x1]
    # corridor: distance to the traced polyline
    yy, xx = np.mgrid[y0:y1, x0:x1]
    d = np.full(sub.shape, np.inf)
    for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
        abx, aby = bx - ax, by - ay
        t = np.clip(((xx - ax) * abx + (yy - ay) * aby) / max(abx * abx + aby * aby, 1e-9), 0, 1)
        d = np.minimum(d, np.hypot(xx - (ax + t * abx), yy - (ay + t * aby)))
    cost = np.where(d <= corridor, 0.02 + sub ** 3, np.inf)
    # vertices to the darkest pixel nearby
    snapped = []
    for x, y in pts.astype(int):
        r0, r1 = max(y - SNAP_PX - y0, 0), y - y0 + SNAP_PX + 1
        c0, c1 = max(x - SNAP_PX - x0, 0), x - x0 + SNAP_PX + 1
        win = sub[r0:r1, c0:c1]
        j, i = np.unravel_index(np.argmin(win), win.shape)
        snapped.append((r0 + j, c0 + i))
    path = []
    for a, b in zip(snapped[:-1], snapped[1:]):
        p, _ = route_through_array(np.where(np.isfinite(cost), cost, 1e9), a, b, fully_connected=True, geometric=True)
        path.extend(p if not path else p[1:])
    return np.array([(x0 + c, y0 + r) for r, c in path], float)       # map pixels (x, y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(DATA.read_text())
    rgb = np.asarray(Image.open(args.map).convert("RGB"))
    v = rgb.max(-1).astype(np.float32) / 255
    n = cs.N
    rid = np.zeros((6, n, n), np.int16)
    seqs = {}
    paths = []
    info = []
    for k, r in enumerate(cfg["rivers"]):
        path = snap_river(v, r["px"], cfg["corridor_px"])
        lat, lon = globe(path[:, 0], path[:, 1])
        f, j, i = cs.latlon_to_cell(lat, lon, n)
        rid[f, j, i] = k + 1
        flat = (f * n * n + j * n + i).astype(np.int64)
        keep = np.r_[True, flat[1:] != flat[:-1]]
        seqs[f"r{k}"] = flat[keep]                 # the river's cells in order, source to mouth
        length = float(np.hypot(*np.diff(path, axis=0).T).sum()) * KM_PER_PX
        paths.append(path)
        info.append({"name": r["name"], "joins": r.get("joins"), "length_km": round(length),
                     "cells": int(len(set(zip(f.tolist(), j.tolist(), i.tolist())))),
                     "source_latlon": [round(math.degrees(float(lat[0])), 3), round(math.degrees(float(lon[0])), 3)],
                     "mouth_latlon": [round(math.degrees(float(lat[-1])), 3), round(math.degrees(float(lon[-1])), 3)]})
    # inland lakes: water pixels (blue, or red for the Bleeding Sea) inside each envelope
    hsv_h, hsv_s, hsv_v = hsv(rgb)
    lid = np.zeros((6, n, n), np.int16)
    for k, lk in enumerate(cfg.get("lakes", [])):
        poly = np.array(lk["polygon_px"])
        x0, y0 = poly.min(0)
        x1, y1 = poly.max(0)
        hh, ss, vv = hsv_h[y0:y1, x0:x1], hsv_s[y0:y1, x0:x1], hsv_v[y0:y1, x0:x1]
        water = (((hh >= 190) & (hh <= 260) & (ss >= 0.06) & (vv >= 0.4))      # blue water
                 | ((ss < 0.10) & (vv > 0.55))                                  # pale grey-blue water (Dagger Lake)
                 | (((hh <= 15) | (hh >= 340)) & (ss >= 0.35)))                 # the red Bleeding Sea
        yy, xx = np.nonzero(water)
        lat, lon = globe(xx + x0 + 0.5, yy + y0 + 0.5)
        f, j, i = cs.latlon_to_cell(lat, lon, n)
        lid[f, j, i] = k + 1
        info.append({"lake": lk["name"], "map_px": int(water.sum()), "cells": int((lid == k + 1).sum()),
                     "area_km2": round(float(water.sum()) * KM_PER_PX ** 2)})
    np.savez_compressed(out / "canon_rivers.npz", river_id=rid, lake_id=lid, **seqs)
    (out / "canon_rivers.json").write_text(json.dumps({"rivers": info}, indent=1, ensure_ascii=False))
    crops = [(800, 1200, 2400, 4100), (2500, 2400, 3800, 4050)]
    panels = []
    for x0, y0, x1, y1 in crops:
        im = Image.fromarray(rgb[y0:y1, x0:x1]).convert("RGB")
        d = ImageDraw.Draw(im)
        for p in paths:
            q = [(x - x0, y - y0) for x, y in p]
            d.line(q, fill=(255, 0, 0), width=3)
        panels.append(im.resize((im.size[0] // 2, im.size[1] // 2)))
    W = sum(p.size[0] for p in panels)
    H = max(p.size[1] for p in panels)
    canvas = Image.new("RGB", (W, H), "white")
    x = 0
    for p in panels:
        canvas.paste(p, (x, 0))
        x += p.size[0]
    canvas.save(out / "check.png")
    print(json.dumps(info, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
