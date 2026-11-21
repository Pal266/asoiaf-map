"""Run the energy balance model (ebm.py) on Earth's land mask and compare with station climate (experiment E09).

Needs the pip package global-land-mask (Earth land/sea at ~1 km, bundled data). No elevation data is used:
station temperatures are reduced to sea level with the lapse rate before comparison.
Usage: python3 earth_validation.py [--out DIR] [--set key=value ...]   (e.g. --set ebm.diffusion_W_m2_K=0.5)
"""

import argparse
import json
from pathlib import Path

import numpy as np
from global_land_mask import globe

import ebm

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"


def earth_land_fraction(res, sub=4):
    lat, lon = ebm.grid(res)
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros((len(lat), len(lon)))
    for a in off:
        for b in off:
            LA, LO = np.meshgrid(np.clip(lat + a * res, -89.99, 89.99), lon + b * res, indexing="ij")
            LO = (LO + 180) % 360 - 180
            acc += globe.is_land(LA, LO)
    return acc / sub ** 2


def sample(field, lat, lon, la, lo):
    """Bilinear sample of a [ny, nx] field (lat descending, periodic lon)."""
    res = lat[0] - lat[1]
    y = (lat[0] - la) / res
    x = (lo - lon[0]) / res
    j0 = int(np.clip(np.floor(y), 0, len(lat) - 2))
    fy = np.clip(y - j0, 0, 1)
    i0 = int(np.floor(x)) % len(lon)
    fx = x - np.floor(x)
    i1 = (i0 + 1) % len(lon)
    return ((1 - fy) * ((1 - fx) * field[j0, i0] + fx * field[j0, i1])
            + fy * ((1 - fx) * field[j0 + 1, i0] + fx * field[j0 + 1, i1]))


def evaluate(cfg, land=None):
    """Run the model on Earth and return (metrics, monthly, ice, lat, lon)."""
    st = json.loads((DATA / "earth_stations.json").read_text())["stations"]
    p = cfg["ebm"]
    land = earth_land_fraction(p["res_deg"]) if land is None else land
    lat, lon, monthly, ice = ebm.run(land, p, cfg["orbit"])
    lapse = cfg["lapse_rate_K_per_km"] / 1000
    A = np.cos(np.radians(lat))[:, None] * np.ones((1, len(lon)))
    rows = []
    for name, la, lo, z, jan, jul, yr in st:
        mj = sample(monthly[0], lat, lon, la, lo) - lapse * z
        ml = sample(monthly[6], lat, lon, la, lo) - lapse * z
        my = sample(monthly.mean(0), lat, lon, la, lo) - lapse * z
        rows.append([name, round(float(mj - jan), 1), round(float(ml - jul), 1), round(float(my - yr), 1)])
    box = json.loads((DATA / "earth_stations.json").read_text()).get("calibration_exclude_box")
    use = np.array([not (box and s[1] > box["lat_min"] and box["lon_min"] < s[2] < box["lon_max"]) for s in st])
    err_all = np.array([r[1:] for r in rows])
    err = err_all[use]
    cellA = (np.cos(np.radians(lat)) * np.radians(lat[0] - lat[1]) * np.radians(lon[1] - lon[0]) * 6371.0 ** 2)[:, None]
    def extent(m, north):
        sel = (lat > 0) if north else (lat < 0)
        return round(float((ice[m] * cellA)[sel].sum()) / 1e6, 2)
    res = {"sea_ice_extent_Mkm2": {"north_mar": extent(2, True), "north_sep": extent(8, True),
                                   "south_sep": extent(8, False), "south_feb": extent(1, False)},
           "sea_ice_extent_earth_Mkm2": {"north_mar": 15.6, "north_sep": 6.5, "south_sep": 18.6, "south_feb": 2.5},
           "stations_used": int(use.sum()), "stations_excluded": [r[0] for r, u in zip(rows, use) if not u],
           "rmse_excluded_jan_jul_year": [round(float(np.sqrt((err_all[~use][:, k] ** 2).mean())), 2) for k in range(3)]
           if (~use).any() else None,
           "global_mean_C": round(float((monthly.mean(0) * A).sum() / A.sum()), 2),
           "rmse_jan_jul_year": [round(float(np.sqrt((err[:, k] ** 2).mean())), 2) for k in range(3)],
           "bias_jan_jul_year": [round(float(err[:, k].mean()), 2) for k in range(3)],
           "zonal_annual_C": {str(round(float(l))): round(float(monthly.mean(0)[j].mean()), 1)
                              for j, l in enumerate(lat) if round(abs(float(l))) in (1, 21, 41, 61, 81, 89)},
           "errors_model_minus_station": rows}
    return res, monthly, ice, lat, lon


def fit_correction(cfg, res, lat, lon, land):
    """Least-squares correction of the model at stations (model statistics, 'MOS'): season-relative residual
    (model - station) = b0 + b1 z_km + b2 M, fitted for the winter and the summer month (January / July swapped in
    the southern hemisphere), with leave-one-out errors."""
    st = json.loads((DATA / "earth_stations.json").read_text())["stations"]
    M = ebm.continentality(lat, lon, land, cfg["ebm"])
    X, Y = [], []
    for (name, la, lo, z, jan, jul, yr), r in zip(st, res["errors_model_minus_station"]):
        if name in res["stations_excluded"]:
            continue
        w, s = (r[1], r[2]) if la >= 0 else (r[2], r[1])
        X.append([1, min(z / 1000, cfg["station_correction"]["z_cap_km"]), sample(M, lat, lon, la, lo)])
        Y.append([w, s])
    X, Y = np.array(X), np.array(Y)
    outc = {}
    for k, lab in enumerate(("winter", "summer")):
        b = np.linalg.lstsq(X, Y[:, k], rcond=None)[0]
        loo = []
        for i in range(len(Y)):
            m = np.ones(len(Y), bool)
            m[i] = False
            bi = np.linalg.lstsq(X[m], Y[m, k], rcond=None)[0]
            loo.append(Y[i, k] - X[i] @ bi)
        outc[lab] = {"coef_1_z_M": [round(float(v), 3) for v in b],
                     "rmse_raw": round(float(np.sqrt((Y[:, k] ** 2).mean())), 2),
                     "rmse_corrected_leave_one_out": round(float(np.sqrt(np.mean(np.square(loo)))), 2)}
    return outc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--set", nargs="*", default=[])
    args = ap.parse_args()
    cfg = json.loads((DATA / "climate.json").read_text())
    for kv in args.set:
        k, v = kv.split("=")
        d = cfg
        parts = k.split(".")
        for q in parts[:-1]:
            d = d[q]
        d[parts[-1]] = float(v)
    land = earth_land_fraction(cfg["ebm"]["res_deg"])
    res, monthly, ice, lat, lon = evaluate(cfg, land)
    res["station_correction_fit"] = fit_correction(cfg, res, lat, lon, land)
    res["params"] = cfg
    print(json.dumps({k: res[k] for k in ("sea_ice_extent_Mkm2", "stations_used", "global_mean_C", "rmse_jan_jul_year", "bias_jan_jul_year",
                                          "rmse_excluded_jan_jul_year", "zonal_annual_C", "station_correction_fit")}))
    worst = sorted(res["errors_model_minus_station"], key=lambda r: -(abs(r[1]) + abs(r[2])))[:12]
    print("worst:", worst)
    if args.out:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "earth_validation.json").write_text(json.dumps(res, indent=1))
        np.savez_compressed(out / "earth_monthly.npz", T=monthly.astype(np.float32), ice=ice, land=land)


if __name__ == "__main__":
    main()
