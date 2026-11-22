"""Test the water balance and routing on Earth (experiment E12).

Earth inputs as in the ice test (world/ice/scripts/earth_ice_validation.py, load_earth): temperature and precipitation
of the E09 and E10 models on Earth, downscaled to 10'; ETOPO1 surface; today's ice. Routing on the 10' grid
(8 neighbours). Each river of data/earth_rivers.json is read at the cell near its gauge (within snap_deg) whose
drainage area is closest to the published basin. Scores: log10 error of discharge, drainage area, and global runoff
to the sea against ~39 000 km3/yr (Dai & Trenberth 2002 give 37 300 km3/yr of continental discharge; with Antarctic
and Greenland ice ~39 000).

Usage: python3 earth_hydro_validation.py --earth-data DIR [--set hydro.fu_w=2.6 ...]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
DATA = HERE.parent / "data"
sys.path.insert(0, str(WORLD / "ice" / "scripts"))
sys.path.insert(0, str(HERE))
import earth_ice_validation as iv     # noqa: E402
import flow                           # noqa: E402

RES = 1 / 6


def latlon_neighbours8(lat, lon):
    ny, nx = len(lat), len(lon)
    idx = np.arange(ny * nx).reshape(ny, nx)
    nbs, dist = [], []
    dy = 111.2e3 * RES
    dx = (np.cos(np.radians(lat)) * 111.2e3 * RES)[:, None] * np.ones((1, nx))
    for dj, di in ((0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)):
        jj = np.clip(np.arange(ny) + dj, 0, ny - 1)
        q = np.roll(idx[jj], -di, 1)
        nbs.append(q.ravel())
        dist.append(np.sqrt((dx * di) ** 2 + (dy * dj) ** 2).ravel() + 1e-9)
    return np.stack(nbs).astype(np.int64), np.stack(dist)


NE_RIVERS = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
             "ne_50m_rivers_lake_centerlines.geojson")


def river_mask(earth_data, lat, lon):
    """Cells crossed by Natural Earth 1:50m river and lake centrelines (public domain), for stream burning: a
    10' grid smoothed by its spherical-harmonic source loses the gorges and valleys that carry rivers (the Congo
    leaves its basin through a narrow gorge), so Earth's big rivers are burnt into it, as in standard DEM
    conditioning (e.g. Saunders 1999). The test is of the water balance, not of Earth's topography."""
    import urllib.request
    f = Path(earth_data) / "ne_50m_rivers_lake_centerlines.geojson"
    if not f.exists():
        urllib.request.urlretrieve(NE_RIVERS, f)
    gj = json.loads(f.read_text())
    mask = np.zeros((len(lat), len(lon)), bool)
    for feat in gj["features"]:
        g = feat["geometry"]
        lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
        for line in lines:
            pts = np.array(line)
            for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
                k = int(max(abs(x1 - x0), abs(y1 - y0)) / (RES / 4)) + 1
                xs = np.linspace(x0, x1, k + 1)
                ys = np.linspace(y0, y1, k + 1)
                j = np.clip(((90 - ys) / RES).astype(int), 0, len(lat) - 1)
                i = (((xs + 180) / RES).astype(int)) % len(lon)
                mask[j, i] = True
    return mask


def evaluate(cfg, earth, nbd=None):
    lat, lon = earth["lat"], earth["lon"]
    ny, nx = len(lat), len(lon)
    nb, dist = nbd if nbd is not None else latlon_neighbours8(lat, lon)
    land = earth["land"].ravel()
    LAT = np.repeat(lat, nx).reshape(ny, nx)
    pet = cfg["hydro"]["pet_factor"] * flow.pet_thornthwaite(earth["T"].astype(np.float64), LAT).sum(0)
    P = earth["P"].sum(0)
    A = ((np.cos(np.radians(lat)) * (111.2 * RES) ** 2)[:, None] * np.ones((1, nx))).ravel()
    z = earth["z"] - cfg["earth_test"]["burn_m"] * earth.get("rivers", 0)
    res = flow.route(z.ravel(), land, P.ravel(), pet.ravel(), A, nb, dist, cfg["hydro"],
                     ice=earth["obs"].ravel(), z_dep=earth["z"].ravel())
    q = res["discharge_m3s"]
    area = flow.accumulate(res["order"], res["rec"], np.where(land, A, 0.0))
    rivers = json.loads((DATA / "earth_rivers.json").read_text())["rivers"]
    rows = []
    LO = np.tile(lon, ny)
    LA = np.repeat(lat, nx)
    for r in rivers:
        near = land & (np.abs(LA - r["lat"]) < cfg["earth_test"]["snap_deg"]) & \
            (np.abs((LO - r["lon"] + 180) % 360 - 180) < cfg["earth_test"]["snap_deg"] / max(np.cos(np.radians(r["lat"])), 0.2))
        cand = np.flatnonzero(near)
        c = cand[np.argmin(np.abs(np.log(np.maximum(area[cand], 1) / r["basin_km2"])))]
        rows.append([r["name"], round(float(q[c])), r["q_m3s"], round(float(area[c])), r["basin_km2"]])
    e = np.array([np.log10(max(x[1], 1) / x[2]) for x in rows])
    ea = np.array([np.log10(max(x[3], 1) / x[4]) for x in rows])
    # runoff reaching the sea: discharge of land cells that drain into sea cells
    rec = res["rec"]
    mouth = land & (rec >= 0) & ~land[np.maximum(rec, 0)]
    to_sea = float(q[mouth].sum()) * flow.SECONDS_PER_YEAR / 1e9
    out = {"rmse_log10_q": round(float(np.sqrt((e ** 2).mean())), 3), "bias_log10_q": round(float(e.mean()), 3),
           "within_factor_1.5": round(float((np.abs(e) < np.log10(1.5)).mean()), 3),
           "rmse_log10_area": round(float(np.sqrt((ea ** 2).mean())), 3),
           "runoff_to_sea_km3_yr": round(to_sea), "terminal_lakes_km2": round(float(A[res["terminal"]].sum())),
           "rivers": rows}
    return out, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--earth-data", required=True)
    ap.add_argument("--set", nargs="*", default=[])
    args = ap.parse_args()
    cfg = json.loads((DATA / "hydrology.json").read_text())
    for kv in args.set:
        k, v = kv.split("=")
        a, b = k.split(".")
        cfg[a][b] = float(v)
    earth = iv.load_earth(args.earth_data)
    earth["rivers"] = river_mask(args.earth_data, earth["lat"], earth["lon"]) & earth["land"]
    out, _ = evaluate(cfg, earth)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
