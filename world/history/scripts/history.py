"""Layer "geological history": crust type, tectonic province, age of the last tectonic event, ocean floor age and
depth on the global cubed sphere (experiment E07).

Input:
  --plates   plates.npz from world/plates/scripts/plates.py (plate_cube, land_cube)
  --relief   canon_relief.npz from world/canon_relief/scripts/canon_relief.py (feature_id)
  world/plates/data/design.json (boundaries) and data/history.json (in the repository)
Output (in --out):
  history.npz   [6, 1000, 1000]: crust (uint8: 0 oceanic, 1 continental), province (uint8, codes in history.json),
                event_age_ma (float32: last tectonic event on continents, floor age in oceans),
                ocean_depth_m (float32: depth of the ocean floor below sea level from its age; NaN on continental crust)
  stats.json    areas, distributions and checks
  preview.png   provinces and ocean floor age, equirectangular 1440 x 720
Deterministic.
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath
from PIL import Image
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
import cubesphere as cs     # noqa: E402

R = cs.R_KM
HISTORY = HERE.parent / "data" / "history.json"
DESIGN = WORLD / "plates" / "data" / "design.json"


def unit(lon_deg, lat_deg):
    lo, la = np.radians(lon_deg), np.radians(lat_deg)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def densify(line, step_deg=0.1):
    pts = []
    for (lo1, la1), (lo2, la2) in zip(line[:-1], line[1:]):
        a, b = unit(lo1, la1), unit(lo2, la2)
        ang = math.acos(max(-1.0, min(1.0, float(a @ b))))
        n = max(2, int(math.degrees(ang) / step_deg) + 1)
        for t in np.linspace(0, 1, n):
            pts.append(a if ang < 1e-9 else (math.sin((1 - t) * ang) * a + math.sin(t * ang) * b) / math.sin(ang))
    return np.array(pts)


def fractal(xyz, spec, rng):
    """Isotropic unit-variance fractal field on the sphere: plane waves along random directions (as in plates.py)."""
    f = np.zeros(xyz.shape[:-1])
    norm = 0.0
    for lam in spec["wavelengths_km"]:
        amp = (lam / spec["wavelengths_km"][0]) ** spec["hurst"]
        for _ in range(spec["components_per_octave"]):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            f += amp * np.cos(2 * np.pi * (R / lam) * (xyz @ d) + rng.uniform(0, 2 * np.pi))
        norm += amp ** 2 * spec["components_per_octave"] / 2
    return f / math.sqrt(norm)


def depth_from_age(t):
    """Basement depth of the ocean floor from its age, plate model GDH1 (Stein and Stein 1992), metres:
    2600 + 365 sqrt(t) below 20 Myr, 5651 - 2473 exp(-0.0278 t) above. No sediment, seamounts or plumes."""
    t = np.maximum(t, 0)
    return np.where(t < 20, 2600 + 365 * np.sqrt(t), 5651 - 2473 * np.exp(-0.0278 * t))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plates", required=True)
    ap.add_argument("--relief", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    hist = json.loads(HISTORY.read_text())
    design = json.loads(DESIGN.read_text())
    bounds = {b["name"]: b for b in design["boundaries"]}

    P = np.load(args.plates)
    plate = P["plate_cube"].astype(np.int16)
    land = P["land_cube"].astype(bool)
    fid = np.load(args.relief)["feature_id"]
    n = cs.N
    lat = np.zeros((6, n, n))
    lon = np.zeros((6, n, n))
    for f in range(6):
        la, lo = cs.cell_latlon(f, n)
        lat[f], lon[f] = np.degrees(la), np.degrees(lo)
    xyz = unit(lon, lat)
    flat = xyz.reshape(-1, 3)

    # ---- distances ----
    cache = {}

    def dist_to_line(name):
        if name not in cache:
            d, _ = cKDTree(densify(bounds[name]["line"])).query(flat)
            cache[name] = (d * R).reshape(lat.shape)
        return cache[name]

    trench_names = [b["name"] for b in design["boundaries"] if b["type"] == "trench"]
    d_trench = np.min([dist_to_line(nm) for nm in trench_names], axis=0)
    land_tree = cKDTree(flat[land.ravel()])
    d_land, nn_land = land_tree.query(flat)
    d_land = (d_land * R).reshape(lat.shape)
    land_index = np.flatnonzero(land.ravel())[nn_land].reshape(lat.shape)

    coast = cs.boundary_cells(land)

    def dist_from(target, query):
        res = np.full(lat.shape, np.inf)
        if target.any() and query.any():
            dd, _ = cKDTree(flat[target.ravel()]).query(flat[query.ravel()])
            res[query] = dd * R
        return res

    # ---- crust ----
    shelf = np.where(d_trench < 300, hist["shelf_width_km"]["active"], hist["shelf_width_km"]["passive"])
    continental = land | (d_land < shelf)
    ocean = ~continental

    # ---- ocean floor age and depth ----
    # each ocean cell takes its age from the nearest ridge that bounds its own plate and could have made it:
    # a ridge reaches at most half_rate x max_age from its axis (floor beyond that is older than the ridge).
    # Plate 0 in a ridge's "between" means every plate it touches (the Antarctic ring). Extinct (fossil) ridges
    # inside a plate date the floor as end_ma + distance / half_rate, up to start_ma. Floor that no ridge reaches
    # takes ocean_default_age_ma.
    ocean_age = np.full(lat.shape, float(hist["ocean_default_age_ma"]))
    best = np.full(lat.shape, np.inf)
    full_rate = np.zeros(lat.shape)          # full spreading rate (mm/yr) of the ridge that made the floor; 0 = none
    for r in hist["ridges"]:
        b = bounds[r["boundary"]]
        d = dist_to_line(r["boundary"])
        plates = [p for p in b["between"] if p]
        on = np.ones(lat.shape, bool) if 0 in b["between"] else np.isin(plate, plates)
        m = on & (d <= r["half_rate_km_per_myr"] * r["max_age_ma"]) & (d < best)
        best[m] = d[m]
        ocean_age[m] = d[m] / r["half_rate_km_per_myr"]
        full_rate[m] = 2 * r["half_rate_km_per_myr"]
    for fr in hist.get("fossil_ridges", []):
        d, _ = cKDTree(densify(fr["line"])).query(flat)
        d = (d * R).reshape(lat.shape)
        m = np.isin(plate, fr["plates"]) & (d <= fr["half_rate_km_per_myr"] * (fr["start_ma"] - fr["end_ma"])) \
            & (d < best)
        best[m] = d[m]
        ocean_age[m] = fr["end_ma"] + d[m] / fr["half_rate_km_per_myr"]
        full_rate[m] = 2 * fr["half_rate_km_per_myr"]
    depth = np.where(ocean, depth_from_age(ocean_age), np.nan).astype(np.float32)

    # oceanic plateaus (large igneous provinces): thick oceanic crust on older floor, away from trenches and land
    op = hist["oceanic_plateaus"]
    prng = np.random.default_rng(hist["seed"] + 1)
    pf = fractal(xyz, op["noise"], prng)
    cand = ocean & (ocean_age > op["min_floor_age_ma"]) & (d_trench > op["min_distance_km"]) \
        & (d_land > op["min_distance_km"])
    Ac = np.broadcast_to(cs.cell_area_km2(n), lat.shape)
    target = op["share_of_ocean"] * Ac[ocean].sum()
    order = np.argsort(-pf[cand])
    cum = np.cumsum(Ac[cand][order])
    thr = pf[cand][order][min(np.searchsorted(cum, target), len(order) - 1)]
    plateau = cand & (pf >= thr)

    # ---- continental provinces ----
    province = np.zeros(lat.shape, np.uint8)
    event_age = np.zeros(lat.shape, np.float32)
    # block outlines are drawn as straight lines; a smooth fractal warp of the lookup coordinates
    # gives them irregular, natural shapes without moving them by more than a few degrees
    rng = np.random.default_rng(hist["seed"])
    amp = hist["warp"]["amplitude_deg"]
    wlon = lon + amp * fractal(xyz, hist["warp"], rng) / np.maximum(np.cos(np.radians(lat)), 0.2)
    wlat = np.clip(lat + amp * fractal(xyz, hist["warp"], rng), -89.99, 89.99)
    wlon = (wlon + 180) % 360 - 180
    pts = np.c_[wlon.ravel(), wlat.ravel()]
    belt_field = fractal(xyz, hist["warp"], rng)
    pts_true = np.c_[lon.ravel(), lat.ravel()]
    for blk in hist["blocks"]:
        # a block with "warp": false is tested in true coordinates (islands that must be hit exactly)
        inside = MPath(np.array(blk["polygon"])).contains_points(
            pts if blk.get("warp", True) else pts_true).reshape(lat.shape)
        m = land & inside & (plate == blk["plate"]) & (province == 0)
        province[m] = blk["province"]
        event_age[m] = blk["age_ma"]
        if "belts" in blk and m.any():
            # cratons are cores welded by younger mobile belts (as Africa's cratons by Pan-African belts)
            thr = np.quantile(belt_field[m], blk["belts"]["share"])
            mb = m & (belt_field <= thr)
            province[mb] = blk["belts"]["province"]
            event_age[mb] = blk["belts"]["age_ma"]
    rest = land & (province == 0)
    province[rest] = 2
    event_age[rest] = 600
    for mg in hist["margins"]:
        # width is counted from the plate's coast that faces the boundary (coast cells within 1000 km of it),
        # not from the boundary line: the schematic E06 lines lie 100-900 km off the generated coasts
        d = dist_to_line(mg["boundary"])
        own = land & (plate == mg["side_plate"])
        facing = own & coast & (d < 1000)
        m = own & (dist_from(facing, own) < mg["width_km"])
        province[m] = mg["province"]
        event_age[m] = mg["age_ma"]
    for k, v in hist["canon_features"].items():
        m = land & (fid == int(k))
        province[m] = v["province"]
        event_age[m] = v["age_ma"]
    # continental shelves take the province of the nearest land
    shelf_cells = continental & ~land
    li = land_index[shelf_cells]
    province[shelf_cells] = province.ravel()[li]
    event_age[shelf_cells] = event_age.ravel()[li]
    event_age[ocean] = ocean_age[ocean]
    province[ocean] = 0
    province[plateau] = 9

    np.savez_compressed(out / "history.npz", crust=continental.astype(np.uint8), province=province,
                        event_age_ma=event_age.astype(np.float32), ocean_depth_m=depth)

    # ---- stats and checks ----
    A = np.broadcast_to(cs.cell_area_km2(n), lat.shape)
    total = float(A.sum())
    names = hist["province_codes"]
    prov_area = {names[str(c)]: round(float(A[province == c].sum())) for c in range(10)}
    oa = ocean_age[ocean]
    ow = A[ocean]
    checks = {
        "oceanic_plateau_share_of_ocean": round(float(A[plateau].sum() / A[ocean].sum()), 4),
        "continental_crust_share": round(float(A[continental].sum()) / total, 4),
        "land_share": round(float(A[land].sum()) / total, 4),
        "ocean_floor_age_ma": {"max": round(float(oa.max()), 1),
                               "area_weighted_mean": round(float((oa * ow).sum() / ow.sum()), 1)},
        "ocean_depth_m": {"area_weighted_mean": round(float((depth[ocean] * ow).sum() / ow.sum())),
                          "min": round(float(depth[ocean].min())), "max": round(float(depth[ocean].max()))},
    }
    # narrow sea floor
    ns = ocean & np.isin(plate, [1, 2]) & (dist_to_line("Narrow sea ridge") < 600) & (lat > 33) & (lat < 60) & (lon < 16)   # between Westeros and Essos
    checks["narrow_sea_floor_max_age_ma"] = round(float(ocean_age[ns].max()), 1) if ns.any() else None
    # oceanic floor by the spreading class of the ridge that made it (Seton et al. 2020 classes, full rate mm/yr)
    oc_floor = ocean & (province != 9)
    classes = {"ultraslow (<20)": (0.1, 20), "slow (20-55)": (20, 55), "intermediate (55-75)": (55, 75),
               "fast (75-180)": (75, 180), "superfast (>180)": (180, 1e9)}
    tot_f = float(A[oc_floor].sum())
    checks["floor_by_spreading_class"] = {k: round(float(A[oc_floor & (full_rate >= a) & (full_rate < b)].sum()) / tot_f, 3)
                                         for k, (a, b) in classes.items()}
    checks["floor_by_spreading_class"]["older than every ridge"] = round(
        float(A[oc_floor & (full_rate == 0)].sum()) / tot_f, 3)
    checks["earth_reference"] = {"floor_by_spreading_class": {"ultraslow": 0.05, "slow": 0.33, "intermediate": 0.20,
                                                              "fast": 0.39, "superfast": 0.03},
                                 "continental_crust_share": 0.41, "mean_ocean_floor_depth_m": 4900,
                                 "max_ocean_floor_age_ma": 180}
    bad = []
    for k, v in hist["canon_features"].items():
        m = land & (fid == int(k))
        if m.any() and not (province[m] == v["province"]).all():
            bad.append(k)
    checks["canon_features_with_wrong_province"] = bad
    stats = {"province_area_km2": prov_area, "checks": checks}
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    # ---- preview ----
    W_, H_ = 1440, 720
    plon = -180 + (np.arange(W_) + 0.5) * 360 / W_
    plat = 90 - (np.arange(H_) + 0.5) * 180 / H_
    LON, LAT = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(LAT), np.radians(LON), n)
    pv = province[ff, jj, ii]
    ag = ocean_age[ff, jj, ii]
    cols = {1: (150, 110, 80), 2: (205, 185, 140), 3: (120, 150, 90), 4: (200, 60, 50), 5: (230, 130, 30),
            6: (180, 90, 200), 7: (230, 210, 90), 8: (90, 170, 160), 9: (120, 200, 230)}
    img = np.zeros((H_, W_, 3), np.uint8)
    oc = pv == 0
    shade = np.clip(1 - ag / 180, 0, 1)
    img[oc] = np.stack([20 + 60 * shade[oc], 40 + 120 * shade[oc], 90 + 140 * shade[oc]], -1).astype(np.uint8)
    for c, col in cols.items():
        img[pv == c] = col
    lnd = land[ff, jj, ii]
    edge = lnd ^ np.roll(lnd, 1, 0) | lnd ^ np.roll(lnd, 1, 1)
    img[edge] = (0, 0, 0)
    Image.fromarray(img).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
