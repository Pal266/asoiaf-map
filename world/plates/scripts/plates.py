"""Present-day plates and the land frame of the whole planet (experiment E06).

Input:
  --outline   canon_outline.npz from world/canon_outline/scripts/canon_outline.py (land on the map)
  data/design.json (in the repository): unknown lands, coast parameters, plate seeds and boundaries
Output (in --out):
  plates.npz       on the tectonic design grid (equirectangular, 0.5 deg, 720 x 360, rows from north to south):
                   land (uint8: 1 land), source (uint8: 0 sea, 1 canon map, 2+k unknown land k), plate (uint8)
                   and on the global cubed sphere [6, 1000, 1000]: plate_cube (uint8), land_cube (uint8)
  plates.geojson   plate boundaries as lines with name and type
  stats.json       areas and checks
  preview.png      plates, land and boundaries (equirectangular, 1440 x 720)
Deterministic: the noise comes from a fixed seed (design.json).
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
import cubesphere as cs     # noqa: E402

R = cs.R_KM
DESIGN = HERE.parent / "data" / "design.json"
RES = 0.5
NLON, NLAT = int(360 / RES), int(180 / RES)


def grid_lonlat():
    lon = -180 + RES * (np.arange(NLON) + 0.5)
    lat = 90 - RES * (np.arange(NLAT) + 0.5)
    return np.meshgrid(lon, lat)


def unit(lon_deg, lat_deg):
    lo, la = np.radians(lon_deg), np.radians(lat_deg)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def cell_area():
    lat_edges = np.radians(90 - RES * np.arange(NLAT + 1))
    band = R * R * np.radians(RES) * (np.sin(lat_edges[:-1]) - np.sin(lat_edges[1:]))
    return np.repeat(band[:, None], NLON, 1)


def densify(line, step_deg=0.25):
    """Points along great circles between successive vertices."""
    pts = []
    for (lo1, la1), (lo2, la2) in zip(line[:-1], line[1:]):
        a, b = unit(lo1, la1), unit(lo2, la2)
        ang = math.acos(max(-1.0, min(1.0, float(a @ b))))
        n = max(2, int(math.degrees(ang) / step_deg) + 1)
        for t in np.linspace(0, 1, n):
            if ang < 1e-9:
                p = a
            else:
                p = (math.sin((1 - t) * ang) * a + math.sin(t * ang) * b) / math.sin(ang)
            pts.append(p)
    return np.array(pts)


def to_lonlat(p):
    return np.degrees(np.arctan2(p[..., 1], p[..., 0])), np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1)))


def inside(poly, lon, lat):
    """Point in polygon in lon/lat (polygons here do not cross the antimeridian)."""
    from matplotlib.path import Path as MPath
    return MPath(np.array(poly)).contains_points(np.c_[lon.ravel(), lat.ravel()]).reshape(lon.shape)


def fractal(xyz, coast, rng):
    """Isotropic fractal field on the sphere: sums of plane waves along random directions (km units)."""
    f = np.zeros(xyz.shape[:-1])
    norm = 0.0
    for lam in coast["wavelengths_km"]:
        amp = (lam / coast["wavelengths_km"][0]) ** coast["hurst"]
        for _ in range(coast["components_per_octave"]):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            ph = rng.uniform(0, 2 * np.pi)
            f += amp * np.cos(2 * np.pi * (R / lam) * (xyz @ d) + ph)
        norm += amp ** 2 * coast["components_per_octave"] / 2
    return f / math.sqrt(norm)


def canon_on_design_grid(outline_path, lon, lat):
    d = np.load(outline_path)
    lf, kf = d["land_frac"], d["known_frac"]
    f, j, i = cs.latlon_to_cell(np.radians(lat), np.radians(lon), cs.N)
    return lf[f, j, i], kf[f, j, i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outline", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    design = json.loads(DESIGN.read_text())
    rng = np.random.default_rng(design["seed"])
    coast = design["coast"]

    lon, lat = grid_lonlat()
    xyz = unit(lon, lat)
    area = cell_area()
    lf, kf = canon_on_design_grid(args.outline, lon, lat)
    known = kf > 0.5
    canon_land = known & (lf > 0.5)

    # ---- unknown lands ----
    noise = fractal(xyz, coast, rng) * coast["amplitude_km"]
    border = known & ~ndimage.binary_erosion(known, iterations=1)
    btree = cKDTree(xyz[border])
    bstate = canon_land[border]
    dist_b, idx_b = btree.query(xyz.reshape(-1, 3))
    dist_b = (dist_b * R).reshape(lon.shape)            # chord ~ arc for short distances
    near_border = ~known & (dist_b < coast["border_blend_km"])
    border_land = bstate[idx_b].reshape(lon.shape)

    # next to the map's edge, unknown cells continue the nearest canon coast
    ext = near_border & border_land
    land = canon_land | ext
    source = np.where(canon_land, 1, 0).astype(np.uint8)
    unknown_stats = []
    for k, u in enumerate(design["unknown_lands"]):
        if "cap_south_of_lat" in u:
            s = (u["cap_south_of_lat"] - lat) * 111.19          # km inside the cap
            zone = lat < u.get("zone_south_of_lat", u["cap_south_of_lat"] + 15)
        else:
            ins = inside(u["polygon"], lon, lat)
            bpts = densify(u["polygon"] + [u["polygon"][0]])
            dist, _ = cKDTree(bpts).query(xyz.reshape(-1, 3))
            dist = (dist * R).reshape(lon.shape)
            s = np.where(ins, dist, -dist)
            zone = s > -1500
        zone &= ~known & (source == 0)
        ext_k = zone & ext
        source[ext_k] = 2 + k                     # the continued coast belongs to this land
        free = zone & ~near_border
        field = s + noise
        need = u["target_km2"] - float(area[ext_k].sum())
        # calibrate the threshold so the land of this frame matches its target area
        lo_t, hi_t = -3000.0, 3000.0
        for _ in range(60):
            t = 0.5 * (lo_t + hi_t)
            if float(area[free & (field > t)].sum()) > need:
                lo_t = t
            else:
                hi_t = t
        m = free & (field > t)
        land |= m
        source[m] = 2 + k
        unknown_stats.append({"name": u["name"], "target_km2": u["target_km2"], "threshold_km": round(t, 1)})
    source[ext & (source == 0)] = 1               # continued coasts outside every frame count with the canon side
    for k, rec in enumerate(unknown_stats):
        rec["area_km2"] = round(float(area[source == 2 + k].sum()))

    # ---- plates ----
    barrier = np.zeros(lon.shape, bool)
    lines_ll = []
    for b in design["boundaries"]:
        p = densify(b["line"], 0.1)
        blo, bla = to_lonlat(p)
        lines_ll.append((b, blo, bla))
        c = np.clip(((blo + 180) / RES).astype(int), 0, NLON - 1)
        r = np.clip(((90 - bla) / RES).astype(int), 0, NLAT - 1)
        barrier[r, c] = True
    barrier = ndimage.binary_dilation(barrier, iterations=1)
    plate = np.zeros(lon.shape, np.uint8)
    open_cells = ~barrier
    # flood fill with wrap-around in longitude: label on a padded array
    pad = np.concatenate([open_cells[:, -2:], open_cells, open_cells[:, :2]], 1)
    lab, _ = ndimage.label(pad)
    lab = lab[:, 2:-2]
    # merge labels connected across the antimeridian
    pairs = set(zip(lab[:, 0][open_cells[:, 0] & open_cells[:, -1]], lab[:, -1][open_cells[:, 0] & open_cells[:, -1]]))
    parent = {}

    def find(a):
        while parent.get(a, a) != a:
            a = parent[a]
        return a
    for a, b in pairs:
        if a and b:
            parent[find(a)] = find(b)
    lab = np.vectorize(lambda a: find(a) if a else 0)(lab)
    seed_label = {}
    for pdef in design["plates"]:
        slo, sla = pdef["seed"]
        c, r = int((slo + 180) / RES), int((90 - sla) / RES)
        seed_label[pdef["id"]] = lab[r, c]
        plate[lab == lab[r, c]] = pdef["id"]
    # cells on barriers or unseeded: nearest assigned cell
    unassigned = plate == 0
    tree = cKDTree(xyz[~unassigned])
    _, nn = tree.query(xyz[unassigned])
    plate[unassigned] = plate[~unassigned][nn]

    # ---- cube grid ----
    n = cs.N
    plate_cube = np.zeros((6, n, n), np.uint8)
    land_cube = np.zeros((6, n, n), np.uint8)
    for face in range(6):
        la, lo = cs.cell_latlon(face, n)
        c = np.clip(((np.degrees(lo) + 180) / RES).astype(int), 0, NLON - 1)
        r = np.clip(((90 - np.degrees(la)) / RES).astype(int), 0, NLAT - 1)
        plate_cube[face] = plate[r, c]
        land_cube[face] = land[r, c]
    np.savez_compressed(out / "plates.npz", land=land.astype(np.uint8), source=source, plate=plate,
                        plate_cube=plate_cube, land_cube=land_cube)

    # ---- stats ----
    total = float(area.sum())
    land_km2 = float(area[land].sum())
    plates_stats = []
    for pdef in design["plates"]:
        m = plate == pdef["id"]
        plates_stats.append({"id": pdef["id"], "name": pdef["name"], "area_km2": round(float(area[m].sum())),
                             "land_km2": round(float(area[m & land].sum()))})
    checks = {}

    def plate_at(lo_, la_):
        return int(plate[int((90 - la_) / RES), int((lo_ + 180) / RES)])
    named = {"Winterfell": (-4.25, 56.51), "King's Landing": (1.18, 38.69), "Sunspear": (4.99, 28.74),
             "Volantis": (22.86, 28.74), "Valyria": (29.7, 22.0), "Qarth": (57.32, 22.46), "Asshai": (83.02, 9.03),
             "Summer Isles (Walano)": (8.0, 17.0), "Ibben": (84.0, 56.0), "Great Moraq": (53.0, 16.0),
             "Dragonmont": (6.16, 42.19), "Basilisk Isles": (34.0, 13.0)}
    checks["plate_of_places"] = {k: plate_at(*v) for k, v in named.items()}
    # distance of Valyria and Dragonmont from the nearest trench / ridge
    trench = np.vstack([unit(blo, bla) for b, blo, bla in lines_ll if "Summer Sea" in b["name"]])
    ridge = np.vstack([unit(blo, bla) for b, blo, bla in lines_ll if b["name"] == "Narrow sea ridge"])
    checks["Valyria_km_from_Summer_trench"] = round(float(cKDTree(trench).query(unit(29.7, 22.0))[0] * R))
    checks["Dragonmont_km_from_narrow_sea_ridge"] = round(float(cKDTree(ridge).query(unit(6.16, 42.19))[0] * R))
    # no land mass may be cut by a ridge or a transform (they lie in oceans); trenches and collisions may border land
    comp, nc = ndimage.label(land)
    split = []
    for c in range(1, nc + 1):
        m = comp == c
        a = float(area[m].sum())
        if a < 1e6:
            continue
        per = {int(p_): round(float(area[m & (plate == p_)].sum())) for p_ in np.unique(plate[m])}
        if len(per) > 1:
            split.append({"area_km2": round(a), "per_plate_km2": per})
    checks["land_masses_on_more_than_one_plate"] = split
    stats = {"land_km2": round(land_km2), "land_share": round(land_km2 / total, 4),
             "canon_land_km2": round(float(area[source == 1].sum())),
             "unknown_lands": unknown_stats, "plates": plates_stats, "checks": checks}
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    gj = [{"type": "Feature", "properties": {"name": b["name"], "type": b["type"], "between": b["between"]},
           "geometry": {"type": "LineString", "coordinates": [[round(float(a), 3), round(float(c), 3)]
                                                              for a, c in zip(blo[::5], bla[::5])]}}
          for b, blo, bla in lines_ll]
    (out / "plates.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": gj}, ensure_ascii=False))

    # ---- preview ----
    pal = {1: (90, 140, 200), 2: (200, 120, 90), 3: (110, 170, 90), 4: (190, 160, 70), 5: (70, 90, 150),
           6: (160, 100, 170), 7: (150, 150, 150)}
    img = np.zeros((NLAT, NLON, 3), np.uint8)
    for pid, col in pal.items():
        img[plate == pid] = col
    img = (img * 0.55).astype(np.uint8)
    img[land] = (img[land] * 0.5 + np.array([235, 220, 170]) * 0.5).astype(np.uint8)
    img[land & (source == 1)] = (img[land & (source == 1)] * 0.35 + np.array([255, 245, 210]) * 0.65).astype(np.uint8)
    prev = Image.fromarray(img).resize((1440, 720), Image.NEAREST)
    dr = ImageDraw.Draw(prev)
    col_t = {"ridge": (255, 255, 0), "trench": (255, 40, 40), "transform": (0, 255, 255)}
    for b, blo, bla in lines_ll:
        xs = (blo + 180) * 4
        ys = (90 - bla) * 4
        for k in range(len(xs) - 1):
            if abs(xs[k + 1] - xs[k]) < 200:
                dr.line([(xs[k], ys[k]), (xs[k + 1], ys[k + 1])], fill=col_t[b["type"]], width=2)
    for pdef in design["plates"]:
        x, y = (pdef["seed"][0] + 180) * 4, (90 - pdef["seed"][1]) * 4
        dr.text((x, y), pdef["name"], fill=(255, 255, 255))
    for gl in range(-60, 61, 30):
        dr.line([(0, (90 - gl) * 4), (1440, (90 - gl) * 4)], fill=(60, 60, 60))
    prev.save(out / "preview.png")
    print(json.dumps({k: v for k, v in stats.items()}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
