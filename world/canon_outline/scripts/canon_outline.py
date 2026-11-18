"""Layer "canon outline": the land/sea outline of the official map on the global cubed-sphere grid.

The process is described in world/canon_outline/GENERATION.md.

Input:  the land mask made by land_mask.py from map M1 (9050 x 6000 px, white = land), and islands_px.json from the
        same run (polygons of the islands smaller than 2 500 km2, in map pixels).
Output (in --out):
  canon_outline.npz  land_frac  float32 [6, 1000, 1000]  share of the cell's known area that is land (0..1)
                     known_frac float32 [6, 1000, 1000]  share of the cell covered by the map (0..1)
  islands.geojson    the islands smaller than 2 500 km2 as polygons in longitude/latitude (degrees), with their
                     area (km2), centroid and the grid cell of the centroid
  validation.json    the checks of GENERATION.md, with their results
  preview.png        equirectangular preview, 3600 x 1800 (sea blue, land beige, outside the map grey)
Deterministic: the same mask and parameters give the same files.
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "grid" / "scripts"))
import cubesphere as cs                                     # noqa: E402
from land_mask import FRAME, LEGEND_BOX, WALL_EAST, WALL_KM, WALL_WEST   # noqa: E402

Image.MAX_IMAGE_PIXELS = None

# Placement of the map on the globe (world/canon/CONSTRAINTS.md, decision 3; experiment E02)
WALL_LAT_DEG = 63.75            # latitude of the Wall's midpoint
WALL_LON_DEG = 0.0              # longitude of the Wall's meridian: a free choice, fixed here
SUB = 8                         # SUB x SUB sample points per cell (~1.2 km apart, about one map pixel)
ROWS_PER_CHUNK = 125

WALL_MID = ((WALL_WEST[0] + WALL_EAST[0]) / 2, (WALL_WEST[1] + WALL_EAST[1]) / 2)
KM_PER_PX = WALL_KM / math.dist(WALL_WEST, WALL_EAST)
R = cs.R_KM
PHI_W = math.radians(WALL_LAT_DEG)
LAM_0 = math.radians(WALL_LON_DEG)


def map_px(lat, lon, turn=0):
    """Sinusoidal projection centred on the Wall: globe (radians) -> map pixel (x, y).
    The map is wider than the planet's circumference north of ~62 N (its north-east corner lies beyond the
    antimeridian), so a point may be found on the map at lon + turn * 360 deg; turn = 0, +1, -1 are tried in turn."""
    lam = (lon - LAM_0 + np.pi) % (2 * np.pi) - np.pi + turn * 2 * np.pi
    x = WALL_MID[0] + R * lam * np.cos(lat) / KM_PER_PX
    y = WALL_MID[1] + R * (PHI_W - lat) / KM_PER_PX
    return x, y


def globe(x, y):
    """Inverse: map pixel -> globe (radians)."""
    lat = PHI_W - (y - WALL_MID[1]) * KM_PER_PX / R
    lon = LAM_0 + (x - WALL_MID[0]) * KM_PER_PX / (R * np.cos(lat))
    return lat, lon


def valid_region(h, w):
    v = np.zeros((h, w), bool)
    v[FRAME:h - FRAME, FRAME:w - FRAME] = True
    x0, y0, x1, y1 = LEGEND_BOX
    v[y0:y1, x0:x1] = False
    return v


def haversine_km(lat1, lon1, lat2, lon2):
    a = math.sin((lat2 - lat1) / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


def face_seam_check(land_frac, known_frac, n):
    """Mean |difference| of land share between cells on a face edge and their neighbours on the next face,
    against the same measure one row inside the face. Similar values mean no seam."""
    across, within = [], []
    e = cs.edges(n)
    d = e[1] - e[0]
    idx = np.arange(n)
    for f in range(6):
        for j, i, dj, di in ((np.zeros(n, int), idx, -1, 0), (np.full(n, n - 1), idx, 1, 0),
                             (idx, np.zeros(n, int), 0, -1), (idx, np.full(n, n - 1), 0, 1)):
            a = e[:-1][i] + (0.5 + di) * d
            b = e[:-1][j] + (0.5 + dj) * d
            lat, lon = cs.xyz_to_latlon(*cs.face_xyz(f, np.tan(a), np.tan(b)))
            g, jj, ii = cs.latlon_to_cell(lat, lon, n)
            m = (known_frac[f, j, i] > 0.99) & (known_frac[g, jj, ii] > 0.99)
            across.append(np.abs(land_frac[f, j, i] - land_frac[g, jj, ii])[m])
            j2 = np.clip(j - dj, 0, n - 1)
            i2 = np.clip(i - di, 0, n - 1)
            m2 = (known_frac[f, j, i] > 0.99) & (known_frac[f, j2, i2] > 0.99)
            within.append(np.abs(land_frac[f, j, i] - land_frac[f, j2, i2])[m2])
    across = np.concatenate(across)
    within = np.concatenate(within)
    return {"pairs": int(across.size), "mean_diff_across_face_edges": round(float(across.mean()), 4),
            "mean_diff_one_row_inside": round(float(within.mean()), 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mask", required=True, help="land_mask.png from land_mask.py")
    ap.add_argument("--islands", required=True, help="islands_px.json from land_mask.py")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    land = np.asarray(Image.open(args.mask)) > 0
    H, W = land.shape
    if (W, H) != (9050, 6000):
        raise SystemExit("expected the 9050 x 6000 land mask from land_mask.py")
    valid = valid_region(H, W)
    labels, _ = ndimage.label(land)
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    top = np.argsort(sizes)[::-1][:6]                        # the six largest land masses on the map
    top_rank = np.zeros(sizes.size, np.int16) - 1
    top_rank[top] = np.arange(top.size)

    n = cs.N
    area = cs.cell_area_km2(n)
    land_frac = np.zeros((6, n, n), np.float32)
    known_frac = np.zeros((6, n, n), np.float32)
    land_area_by_mass = np.zeros(top.size)
    e = cs.edges(n)
    d = e[1] - e[0]
    offs = (np.arange(SUB) + 0.5) / SUB

    for f in range(6):
        for j0 in range(0, n, ROWS_PER_CHUNK):
            j1 = min(n, j0 + ROWS_PER_CHUNK)
            lsum = np.zeros((j1 - j0, n))
            ksum = np.zeros((j1 - j0, n))
            for oj in offs:
                b = e[j0:j1] + oj * d
                for oi in offs:
                    a = e[:-1] + oi * d
                    t, u = np.meshgrid(np.tan(a), np.tan(b))
                    lat, lon = cs.xyz_to_latlon(*cs.face_xyz(f, t, u))
                    xi = np.zeros(lat.shape, np.int64)
                    yi = np.zeros(lat.shape, np.int64)
                    k = np.zeros(lat.shape, bool)
                    for turn in (0, 1, -1):
                        x, y = map_px(lat, lon, turn)
                        xt = np.floor(x).astype(np.int64)
                        yt = np.floor(y).astype(np.int64)
                        ok = ~k & (xt >= 0) & (xt < W) & (yt >= 0) & (yt < H)
                        ok[ok] = valid[yt[ok], xt[ok]]
                        xi[ok] = xt[ok]
                        yi[ok] = yt[ok]
                        k |= ok
                    l = np.zeros(xi.shape, bool)
                    l[k] = land[yi[k], xi[k]]
                    ksum += k
                    lsum += l
                    rank = np.full(xi.shape, -1, np.int16)
                    rank[l] = top_rank[labels[yi[l], xi[l]]]
                    a_cell = area[j0:j1, :] / SUB ** 2
                    for r in range(top.size):
                        land_area_by_mass[r] += a_cell[rank == r].sum()
            known_frac[f, j0:j1] = ksum / SUB ** 2
            land_frac[f, j0:j1] = np.where(ksum > 0, lsum / np.maximum(ksum, 1), 0)

    np.savez_compressed(out / "canon_outline.npz", land_frac=land_frac, known_frac=known_frac)

    # ---- island polygons in longitude / latitude ----
    islands_px = json.loads(Path(args.islands).read_text())
    features = []
    island_grid_km2 = 0.0
    for isl in islands_px:
        px = np.array(isl["polygon_px"], float)
        la, lo = globe(px[:, 0] + 0.5, px[:, 1] + 0.5)          # pixel centres
        lo = (lo + np.pi) % (2 * np.pi) - np.pi
        cla, clo = globe(isl["centroid_px"][0] + 0.5, isl["centroid_px"][1] + 0.5)
        clo = (clo + np.pi) % (2 * np.pi) - np.pi
        f_, j_, i_ = cs.latlon_to_cell(np.array([cla]), np.array([clo]), n)
        ring = [[round(math.degrees(b), 5), round(math.degrees(a), 5)] for a, b in zip(la, lo)]
        ring.append(ring[0])
        features.append({"type": "Feature",
                         "properties": {"id": isl["id"], "area_km2": isl["area_km2"],
                                        "centroid_lon_lat": [round(math.degrees(clo), 4), round(math.degrees(cla), 4)],
                                        "cell_face_j_i": [int(f_[0]), int(j_[0]), int(i_[0])]},
                         "geometry": {"type": "Polygon", "coordinates": [ring]}})
        island_grid_km2 += isl["area_km2"]
    (out / "islands.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": features}))

    # ---- validation ----
    A = np.broadcast_to(area, (6, n, n))
    grid_known_km2 = float((known_frac * A).sum())
    grid_land_km2 = float((land_frac * known_frac * A).sum())
    # the same quantities measured directly on the map, which is equal-area under this projection
    px_km2 = KM_PER_PX ** 2
    map_known_km2 = float(valid.sum() * px_km2)
    map_land_km2 = float((land & valid).sum() * px_km2)
    mass_map_km2 = sizes[top] * px_km2

    lw = globe(*WALL_WEST)
    le = globe(*WALL_EAST)
    wall_km = haversine_km(lw[0], lw[1], le[0], le[1])

    named = {"Winterfell": (1550, 1607), "Castle Black": (1798, 987), "Sunspear": (2124, 3980),   # map px, read on M1
             "Oldtown": (930, 3873), "King's Landing": (1829, 3130), "Volantis": (3463, 3980),
             "Qarth": (6276, 4516), "Asshai": (8756, 5664)}
    places = {}
    for name, (px, py) in named.items():
        la, lo = globe(px, py)
        places[name] = {"lat": round(math.degrees(la), 2), "lon": round(math.degrees(lo), 2)}

    total_area = float(A.sum())
    checks = {
        "total_grid_area_km2": round(total_area),
        "sphere_area_km2": round(4 * math.pi * R ** 2),
        "known_area_km2": {"grid": round(grid_known_km2), "map": round(map_known_km2),
                           "rel_diff": round(abs(grid_known_km2 - map_known_km2) / map_known_km2, 5)},
        "land_area_km2": {"grid": round(grid_land_km2), "map": round(map_land_km2),
                          "rel_diff": round(abs(grid_land_km2 - map_land_km2) / map_land_km2, 5)},
        "largest_land_masses_km2": [
            {"map": round(float(m)), "grid": round(float(g)), "rel_diff": round(abs(g - m) / m, 4)}
            for m, g in zip(mass_map_km2, land_area_by_mass)],
        "wall_length_km": round(wall_km, 1),
        "islands_below_2500_km2": {"count": len(features), "area_km2": round(island_grid_km2, 1),
                                   "smaller_than_one_cell_100_km2": sum(1 for f in features
                                                                        if f["properties"]["area_km2"] < 100)},
        "face_seams": face_seam_check(land_frac, known_frac, n),
        "known_share_of_planet": round(grid_known_km2 / total_area, 4),
        "land_share_of_known": round(grid_land_km2 / grid_known_km2, 4),
        "places": places,
        "parameters": {"wall_lat_deg": WALL_LAT_DEG, "wall_lon_deg": WALL_LON_DEG,
                       "km_per_px": round(KM_PER_PX, 4), "sub_samples": SUB ** 2},
    }
    (out / "validation.json").write_text(json.dumps(checks, indent=2, ensure_ascii=False))

    # ---- preview ----
    pw, ph = 3600, 1800
    lon = np.radians((np.arange(pw) + 0.5) / pw * 360 - 180)
    lat = np.radians(90 - (np.arange(ph) + 0.5) / ph * 180)
    LON, LAT = np.meshgrid(lon, lat)
    ff, jj, ii = cs.latlon_to_cell(LAT, LON, n)
    kf = known_frac[ff, jj, ii]
    lf = land_frac[ff, jj, ii]
    img = np.full((ph, pw, 3), 128, np.uint8)
    sea = np.array([40, 70, 130], float)
    beige = np.array([225, 205, 160], float)
    col = sea * (1 - lf[..., None]) + beige * lf[..., None]
    grey = np.array([128, 128, 128], float)
    img = (grey * (1 - kf[..., None]) + col * kf[..., None]).astype(np.uint8)
    for gl in range(-60, 61, 30):                            # graticule
        r = int((90 - gl) / 180 * ph)
        img[r, :] = (img[r, :] * 0.6).astype(np.uint8)
    for gm in range(-180, 180, 30):
        c = int((gm + 180) / 360 * pw)
        img[:, c] = (img[:, c] * 0.6).astype(np.uint8)
    Image.fromarray(img).save(out / "preview.png")
    print(json.dumps(checks, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
