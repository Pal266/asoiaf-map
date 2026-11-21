"""Layer "relief": mean elevation and sub-grid relief of every cell of the global cubed sphere (experiment E08).

Input:
  --plates   plates.npz (plate_cube, land_cube)                   from world/plates/scripts/plates.py
  --relief   canon_relief.npz (feature_id, relief)                from world/canon_relief/scripts/canon_relief.py
  --history  history.npz (crust, province, ocean_depth_m)         from world/history/scripts/history.py
  world/plates/data/design.json and world/relief/data/relief.json (in the repository)
Output (in --out):
  relief.npz   [6, 1000, 1000] float32: elevation_m (cell mean, metres above sea level, negative = sea floor),
               subgrid_relief_m (expected peak-to-valley relief inside the cell)
  stats.json   hypsometry and checks
  preview.png  shaded elevation, equirectangular 1440 x 720
Deterministic. Bedrock surface: no ice, no lakes, no river incision (later layers).
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
import cubesphere as cs     # noqa: E402

R = cs.R_KM
PARAMS = HERE.parent / "data" / "relief.json"
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
    """Isotropic unit-variance fractal field on the sphere (plane waves along random directions)."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plates", required=True)
    ap.add_argument("--relief", required=True)
    ap.add_argument("--history", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    P = json.loads(PARAMS.read_text())
    design = json.loads(DESIGN.read_text())
    bounds = {b["name"]: b for b in design["boundaries"]}

    pl = np.load(args.plates)
    plate = pl["plate_cube"].astype(np.int16)
    land = pl["land_cube"].astype(bool)
    cr = np.load(args.relief)
    fid, dens = cr["feature_id"], cr["relief"]
    hs = np.load(args.history)
    crust, province, odepth = hs["crust"].astype(bool), hs["province"], hs["ocean_depth_m"]

    n = cs.N
    lat = np.zeros((6, n, n))
    lon = np.zeros((6, n, n))
    for f in range(6):
        la, lo = cs.cell_latlon(f, n)
        lat[f], lon[f] = np.degrees(la), np.degrees(lo)
    xyz = unit(lon, lat)
    flat = xyz.reshape(-1, 3)
    shape = lat.shape

    def dist_to_line(name):
        d, _ = cKDTree(densify(bounds[name]["line"])).query(flat)
        return (d * R).reshape(shape)

    def dist_to_cells(target, query_mask):
        """Distance (km) from cells of query_mask to the nearest cell of target; inf elsewhere."""
        res = np.full(shape, np.inf)
        if target.any() and query_mask.any():
            d, idx = cKDTree(flat[target.ravel()]).query(flat[query_mask.ravel()])
            res[query_mask] = d * R
        return res

    def interior(mask, half_width_km=None):
        """Normalised depth inside a region: 0 at its edge, 1 at the crest (max depth or half_width)."""
        d = dist_to_cells(cs.boundary_cells(mask, halo), mask)
        d[~mask] = 0
        scale = half_width_km if half_width_km else max(float(d[mask].max()), 1.0)
        return np.clip(d / scale, 0, 1)

    halo = cs.halo_index(n, 1)
    rng = np.random.default_rng(P["seed"])
    noise_plain = fractal(xyz, P["noise"], rng)
    noise_range = fractal(xyz, P["noise"], rng)
    noise_sea = fractal(xyz, P["noise"], rng)
    ridged = np.clip(1 - np.abs(fractal(xyz, P["range_noise"], rng)), 0, 1)
    residual = P["ocean"]["residual_depth"]["rms_m"] * fractal(xyz, P["ocean"]["residual_depth"], rng)

    trench_names = [b["name"] for b in design["boundaries"] if b["type"] == "trench"]
    d_trench = np.min([dist_to_line(nm) for nm in trench_names], axis=0)
    ocean = ~land
    d_coast = dist_to_cells(ocean, land)            # land cells: distance to the sea
    d_land = dist_to_cells(land, ocean)             # sea cells: distance to land
    active = d_trench < 300

    # ---------------- land ----------------
    mean = np.zeros(shape)
    amp = np.zeros(shape)
    for code, v in P["plain_m"].items():
        m = province == int(code)
        mean[m], amp[m] = v["mean"], v["amplitude"]
    # province levels change gradually across block borders (sigma ~ 50 km), not in steps
    mean = cs.box_smooth(mean, P["plain_smoothing_iterations"], halo)
    amp = cs.box_smooth(amp, P["plain_smoothing_iterations"], halo)
    plain = np.maximum(mean + amp * noise_plain, 0.25 * mean)
    W = np.where(active, P["coastal_plain_km"]["active"], P["coastal_plain_km"]["passive"])
    ramp = np.clip(np.where(land, d_coast, 0) / W, 0, 1) ** 0.7
    elev = 2 + (plain - 2) * ramp

    crest = np.zeros(shape)       # metres above the plain at the crest
    prof = np.zeros(shape)        # 0 at the edge of a range, 1 on its crest
    cf = P["crest_factor"]
    canon_range_ids = [int(k) for k in P["canon_ranges"]]
    generic_ok = land & ~np.isin(fid, canon_range_ids)
    for code, g in P["generic_ranges"].items():
        m = generic_ok & (province == int(code))
        if m.any():
            crest[m] = cf * g["peak_m"]
            prof[m] = interior(m, g["half_width_km"])[m]
    for mg in P["margin_ranges"]:
        d = dist_to_line(mg["boundary"])
        own = land & (plate == mg["side_plate"])
        facing = own & cs.boundary_cells(land, halo) & (d < 1000)
        m = generic_ok & own & (dist_to_cells(facing, own) < mg["width_km"]) & np.isin(province, [4, 5])
        if m.any():
            crest[m] = cf * mg["peak_m"]
            prof[m] = interior(m, mg["half_width_km"])[m]
    for k, c in P["canon_ranges"].items():
        m = land & (fid == int(k))
        if m.any():
            drawn = np.clip(dens[m] / 0.5, 0, 1)          # drawn relief density of the official map
            crest[m] = cf * c["peak_m"] * (0.7 + 0.3 * drawn)
            prof[m] = interior(m)[m]
    peak_term = np.zeros(shape)
    for k, c in P["canon_peaks"].items():
        if k == "note":
            continue
        d = np.linalg.norm(xyz - unit(c["lon"], c["lat"]), axis=-1) * R
        if c.get("shape") == "island":           # a ring island: its own land takes a range profile
            m = land & (d < c["radius_km"])
            cone = np.zeros(shape)
            cone[m] = cf * c["peak_m"] * np.sin(0.5 * np.pi * interior(m)[m])
            peak_term = np.maximum(peak_term, cone)
            continue
        cone = cf * c["peak_m"] * np.clip(1 - d / c["radius_km"], 0, 1)
        if not (land & (cone > 0)).any():          # a peak narrower than a cell still lifts its cell
            f, j, i = cs.latlon_to_cell(np.radians(c["lat"]), np.radians(c["lon"]), n)
            cone[f, j, i] = 0.3 * cf * c["peak_m"]
        peak_term = np.maximum(peak_term, np.where(land, cone, 0))
    rng_term = crest * np.sin(0.5 * np.pi * prof) ** P["profile_exponent"] * np.clip(0.8 + 0.25 * noise_range, 0.4, 1.3) \
        * (0.5 + 0.8 * ridged)
    # round the hand-drawn envelope corners (sigma ~ 26 km)
    rng_term = cs.box_smooth(rng_term, P["range_smoothing_iterations"], halo)
    rng_term *= np.clip(np.where(land, d_coast, 0) / 10.0, 0, 1)
    rng_term = np.maximum(rng_term, peak_term)
    elev = elev + rng_term
    basin_term = np.zeros(shape)
    for k, b in P["canon_basins"].items():
        m = land & (fid == int(k))
        if m.any():
            basin_term[m] = b["depth_m"] * np.sin(0.5 * np.pi * interior(m)[m])
    basin_term = cs.box_smooth(basin_term, P["range_smoothing_iterations"], halo)
    elev = elev - basin_term
    for k, lw in P["canon_lowlands"].items():
        m = land & (fid == int(k))
        w = np.zeros(shape)
        w[m] = np.sin(0.5 * np.pi * interior(m)[m])
        w = cs.box_smooth(w, P["range_smoothing_iterations"], halo)
        elev = elev - w * np.maximum(elev - lw["max_m"], 0)
    elev = np.where(land, np.maximum(elev, 1.0), elev)

    # ---------------- sea floor ----------------
    oc = P["ocean"]
    oceanic = ocean & ~crust
    margin = ocean & crust
    floor = -odepth.astype(np.float64)
    sd = oc["sediment"]
    ages = (odepth.astype(np.float64) - 2600).clip(0)      # invert GDH1 only to recover age for sediment
    age = np.where(odepth < 2600 + 365 * np.sqrt(20), (ages / 365) ** 2,
                   -np.log(np.clip((5651 - odepth.astype(np.float64)) / 2473, 1e-6, None)) / 0.0278)
    rho = sd["density_kg_m3"]
    sed = np.minimum(sd["pelagic_m_per_myr"] * np.nan_to_num(age), sd["pelagic_max_m"]) \
        + sd["apron_m"] * np.exp(-(np.where(ocean, d_land, 0) / sd["apron_scale_km"]) ** 2)
    floor = floor + sed * (rho["mantle"] - rho["sediment"]) / (rho["mantle"] - rho["water"])
    tr = oc["trench_extra_m"] * np.exp(-(d_trench / oc["trench_halfwidth_km"]) ** 2)
    floor = floor - tr + oc["abyssal_hill_amplitude_m"] * noise_sea
    # oceanic plateaus stand ~2 km above the floor (Ontong Java, Kerguelen), with blurred flanks
    plat = cs.box_smooth((province == 9).astype(float), P["range_smoothing_iterations"], halo)
    floor = floor + oc["plateau_height_m"] * plat + residual
    # continental margin: shelf, then slope down to the depth of the nearest oceanic floor
    tree_oc = cKDTree(flat[oceanic.ravel()])
    _, ix = tree_oc.query(flat[margin.ravel()])
    abyss = np.zeros(shape)
    abyss[margin] = -floor[oceanic][ix]
    sw = json.loads((WORLD / "history" / "data" / "history.json").read_text())["shelf_width_km"]
    Wm = np.where(active, float(sw["active"]), float(sw["passive"]))
    shelf_w = oc["shelf_share_of_margin"] * Wm
    dl = np.where(margin, d_land, 0)
    t = np.clip((dl - shelf_w) / (Wm - shelf_w), 0, 1)
    smooth = t * t * (3 - 2 * t)
    shelf_depth = np.where(dl < shelf_w, 5 + (oc["shelf_break_m"] - 5) * dl / shelf_w,
                           oc["shelf_break_m"] + (abyss - oc["shelf_break_m"]) * smooth)
    sea = np.where(margin, -shelf_depth, floor)
    for k, s in P["canon_sunken"].items():
        m = ocean & (fid == int(k))
        sea[m] = -s["depth_m"]
    elev = np.where(land, elev, np.minimum(sea, -1.0))

    sg = P["subgrid_relief"]
    subgrid = np.where(land, sg["base_m"] + sg["per_metre_of_range"] * rng_term + 0.5 * basin_term,
                       2 * oc["abyssal_hill_amplitude_m"])
    sediment = np.where(ocean, sed, 0)
    np.savez_compressed(out / "relief.npz", elevation_m=elev.astype(np.float32),
                        subgrid_relief_m=subgrid.astype(np.float32), sediment_m=sediment.astype(np.float32))

    # ---------------- stats ----------------
    A = np.broadcast_to(cs.cell_area_km2(n), shape)
    la_, oa_ = A[land], A[ocean]
    el, eo = elev[land], elev[ocean]

    def wq(v, w, q):
        o = np.argsort(v)
        c = np.cumsum(w[o]) / w.sum()
        return float(v[o][np.searchsorted(c, q)])

    feats = {}
    for k, c in list(P["canon_ranges"].items()) + [(k, v) for k, v in P["canon_peaks"].items() if k != "note"]:
        m = land & (fid == int(k))
        if k in P["canon_peaks"]:
            d = np.linalg.norm(xyz - unit(c["lon"], c["lat"]), axis=-1) * R
            m = land & (d <= max(c["radius_km"], 10))
        if m.any():
            feats[k] = {"analogue_peak_m": c["peak_m"], "highest_cell_mean_m": round(float(elev[m].max())),
                        "feature_mean_m": round(float((elev[m] * A[m]).sum() / A[m].sum()))}
    # sediment budget against GlobSed: margin zone = the 12.9 % of the ocean nearest land
    dlo = d_land[ocean]
    o_order = np.argsort(dlo)
    cumA = np.cumsum(oa_[o_order]) / oa_.sum()
    Dm = float(dlo[o_order][np.searchsorted(cumA, 0.129)])
    mz, far = ocean & (d_land <= Dm), ocean & (d_land > Dm + 200)
    sed_budget = {
        "margin_zone_width_km": round(Dm),
        "margin_zone_mean_m": round(float((sediment[mz] * A[mz]).sum() / A[mz].sum())),
        "far_ocean_mean_m": round(float((sediment[far] * A[far]).sum() / A[far].sum())),
        "total_million_km3": round(float((sediment[ocean] * A[ocean]).sum()) / 1e9, 1),
        "globsed": {"margin_zone_mean_m": 3044, "far_ocean_mean_m": 404, "total_million_km3": 337}}
    stats = {
        "sediment": sed_budget,
        "land_mean_elevation_m": round(float((el * la_).sum() / la_.sum())),
        "land_elevation_quantiles_m": {q: round(wq(el, la_, q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99)},
        "land_share_above_m": {h: round(float(la_[el > h].sum() / la_.sum()), 4) for h in (500, 1000, 2000, 3000)},
        "ocean_mean_depth_m": round(float(-(eo * oa_).sum() / oa_.sum())),
        "ocean_depth_quantiles_m": {q: round(-wq(eo, oa_, q)) for q in (0.1, 0.5, 0.9)},
        "highest_cell_m": round(float(el.max())), "deepest_cell_m": round(float(eo.min())),
        "coast_consistent": bool(((elev > 0) == land).all()),
        "earth_reference": {"land_mean_elevation_m": 840, "ocean_mean_depth_m": 3688},
        "canon_ranges": feats,
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    # ---------------- preview ----------------
    W_, H_ = 1440, 720
    plon = -180 + (np.arange(W_) + 0.5) * 360 / W_
    plat = 90 - (np.arange(H_) + 0.5) * 180 / H_
    LON, LAT = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(LAT), np.radians(LON), n)
    e = elev[ff, jj, ii]
    gy, gx = np.gradient(e)
    shade = np.clip(1 - 0.004 * (gx - gy), 0.6, 1.3)
    img = np.zeros((H_, W_, 3))
    s = np.clip(-e / 6000, 0, 1)
    img[e < 0] = np.stack([30 + 60 * (1 - s), 60 + 110 * (1 - s), 120 + 110 * (1 - s)], -1)[e < 0]
    stops = np.array([[0, 70, 140, 60], [300, 140, 180, 90], [800, 210, 200, 120], [1800, 170, 120, 70],
                      [3500, 230, 230, 230], [6000, 255, 255, 255]], float)
    for c in range(3):
        img[..., c] = np.where(e >= 0, np.interp(e, stops[:, 0], stops[:, c + 1]), img[..., c])
    img[e >= 0] *= shade[e >= 0, None]
    Image.fromarray(np.clip(img, 0, 255).astype(np.uint8)).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
