"""Run the moisture model (moisture.py) on Earth and compare with station precipitation (experiment E10).

Earth inputs: land mask (pip package global-land-mask), topography at 10 arc-minutes (ETOPO1 ice surface, from
fatiando-data/earth-topography-10arcmin, DOI 10.5281/zenodo.5882203; downloaded once into --earth-data), and
monthly surface temperature from the temperature model (world/climate/scripts/ebm.py on Earth).
Usage: python3 earth_precip_validation.py --earth-data DIR [--out DIR] [--set key=value ...]
"""

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
DATA = HERE.parent / "data"
sys.path.insert(0, str(WORLD / "climate" / "scripts"))
import ebm                                    # noqa: E402
import earth_validation as tval               # noqa: E402  (the temperature test: land mask, sampling)
import moisture                               # noqa: E402

TOPO_URL = "https://github.com/fatiando-data/earth-topography-10arcmin/releases/download/v1/earth-topography-10arcmin.nc"


def earth_topography(path, res):
    """Mean elevation of land per res-degree cell (m; 0 at sea), from the 10' grid."""
    import h5py
    f = path / "earth-topography-10arcmin.nc"
    if not f.exists():
        path.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(TOPO_URL, f)
    h = h5py.File(f, "r")
    z = h["topography"][:].astype(float)[:-1, :-1]            # drop the duplicate last row/column of nodes
    z = z[::-1]                                                  # north first
    k = int(round(res * 6))
    ny, nx = z.shape[0] // k, z.shape[1] // k
    zb = z[:ny * k, :nx * k].reshape(ny, k, nx, k)
    return np.where(zb > 0, zb, 0).mean((1, 3))


def to_grid(field, lat_c, lon_c, lat, lon):
    """Bilinear interpolation from a coarse (lat_c, lon_c) grid to (lat, lon)."""
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    out = np.zeros(LA.shape)
    for j in range(len(lat)):
        for i in range(len(lon)):
            out[j, i] = tval.sample(field, lat_c, lon_c, LA[j, i], LO[j, i])
    return out


def land_sample(field, land, lat, lon, la, lo):
    """Bilinear interpolation weighted by land fraction, so a coastal station reads the land cells
    around it rather than the open ocean next door (stations are on land; the ocean cells of a coarse
    grid carry the marine rain rate). Falls back to plain bilinear when no neighbour has land."""
    dl = lat[1] - lat[0]
    fi = (la - lat[0]) / dl
    fj = ((lo - lon[0]) / (lon[1] - lon[0])) % len(lon)
    i0 = int(np.clip(np.floor(fi), 0, len(lat) - 2)); j0 = int(np.floor(fj)) % len(lon)
    ti = np.clip(fi - i0, 0, 1); tj = fj - np.floor(fj)
    num = den = num0 = 0.0
    for di, wi in ((0, 1 - ti), (1, ti)):
        for dj, wj in ((0, 1 - tj), (1, tj)):
            ii, jj = i0 + di, (j0 + dj) % len(lon)
            w = wi * wj
            num0 += w * field[ii, jj]
            num += w * land[ii, jj] * field[ii, jj]; den += w * land[ii, jj]
    return num / den if den > 1e-6 else num0


def evaluate(cfg, earth):
    p = cfg
    lat, lon = moisture.grid(p["precip"]["res_deg"])
    P, E, w, diag = moisture.run(lat, lon, earth["T"], earth["land"], earth["elev"], earth["ice"], p)
    st = json.loads((DATA / "earth_precip_stations.json").read_text())["stations"]
    rows = []
    for name, la, lo, yr, jan, jul in st:
        my = land_sample(P.sum(0), earth["land"], lat, lon, la, lo)
        mj = land_sample(P[0], earth["land"], lat, lon, la, lo)
        ml = land_sample(P[6], earth["land"], lat, lon, la, lo)
        rows.append([name, round(float(my)), round(float(mj)), round(float(ml)), yr, jan, jul])
    r = np.array([x[1:] for x in rows], float)
    lg = np.log10(np.maximum(r[:, 0], 5)) - np.log10(np.maximum(r[:, 3], 5))
    share_m = (r[:, 1] - r[:, 2]) / np.maximum(r[:, 0], 1)        # (Jan - Jul) / year: the seasonal contrast
    share_o = (r[:, 4] - r[:, 5]) / np.maximum(r[:, 3], 1)
    A = np.cos(np.radians(lat))[:, None] * np.ones((1, len(lon)))
    res = {"global_mean_mm": round(float((P.sum(0) * A).sum() / A.sum())),
           "land_mean_mm": round(float((P.sum(0) * A * earth["land"]).sum() / (A * earth["land"]).sum())),
           "rmse_log10_annual": round(float(np.sqrt((lg ** 2).mean())), 3),
           "bias_log10_annual": round(float(lg.mean()), 3),
           "rmse_seasonal_contrast": round(float(np.sqrt(((share_m - share_o) ** 2).mean())), 3),
           "within_factor_2": round(float((np.abs(lg) < np.log10(2)).mean()), 3),
           "stations": rows}
    obs = r[:, 3]
    arid = obs < 250
    res["arid_bias_log10"] = round(float(lg[arid].mean()), 3)                 # stations below 250 mm a year
    zon = json.loads((DATA / "earth_zonal_gpcp.json").read_text())
    Z = P.sum(0).mean(1)
    zm = np.array([Z[np.argmin(np.abs(lat - la))] for la in zon["lat_deg"]])
    res["zonal_mm"] = [int(round(z)) for z in zm]
    res["zonal_rmse_log10"] = round(float(np.sqrt(np.mean((np.log10(zm) - np.log10(zon["mm_per_year"])) ** 2))), 3)
    res["_diag"] = diag
    return res, P, E


def load_earth(cfg, earth_data):
    lat, lon = moisture.grid(cfg["precip"]["res_deg"])
    tcfg = json.loads((WORLD / "climate" / "data" / "climate.json").read_text())
    land2 = tval.earth_land_fraction(tcfg["ebm"]["res_deg"])
    lat2, lon2, Tm, ice = ebm.run(land2, tcfg["ebm"], tcfg["orbit"])
    T = np.stack([to_grid(Tm[m], lat2, lon2, lat, lon) for m in range(12)])
    icef = np.stack([to_grid(ice[m].astype(float), lat2, lon2, lat, lon) > 0.5 for m in range(12)])
    land = tval.earth_land_fraction(cfg["precip"]["res_deg"])
    elev = earth_topography(Path(earth_data), cfg["precip"]["res_deg"])
    return {"T": T, "ice": icef, "land": land, "elev": elev}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--earth-data", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--set", nargs="*", default=[])
    args = ap.parse_args()
    cfg = json.loads((DATA / "precipitation.json").read_text())
    for kv in args.set:
        k, v = kv.split("=")
        d = cfg
        parts = k.split(".")
        for q in parts[:-1]:
            d = d[q]
        d[parts[-1]] = float(v)
    earth = load_earth(cfg, args.earth_data)
    res, P, E = evaluate(cfg, earth)
    res.pop("_diag")
    print(json.dumps({k: v for k, v in res.items() if k != "stations"}))
    worst = sorted(res["stations"], key=lambda r: -abs(np.log10(max(r[1], 5)) - np.log10(max(r[4], 5))))[:15]
    print("worst [name, model year, jan, jul, station year, jan, jul]:", worst)
    if args.out:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "earth_validation.json").write_text(json.dumps(res, indent=1))
        np.savez_compressed(out / "earth_precip.npz", P=P.astype(np.float32), **{k: v for k, v in earth.items()})


if __name__ == "__main__":
    main()
