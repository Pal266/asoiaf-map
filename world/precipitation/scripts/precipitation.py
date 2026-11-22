"""Layer "climate: precipitation": monthly precipitation of every cell of the global cubed sphere (experiment E10).

The moisture-budget model (moisture.py, calibrated on Earth by earth_precip_validation.py) runs on a 2-degree grid (the resolution it was calibrated at)
driven by the temperature model (world/climate/scripts/ebm.py with world/climate/data/climate.json, the same
run as the temperature layer) and the world's relief; monthly totals are interpolated to the 10 km grid.

Input:
  --relief          relief.npz (elevation_m)           from world/relief/scripts/relief.py
  --canon-relief    canon_relief.npz (feature_id)       from world/canon_relief/scripts/canon_relief.py
  world/precipitation/data/precipitation.json, world/climate/data/climate.json (in the repository)
Output (in --out):
  precipitation.npz  p_month_mm uint16 [12, 6, 1000, 1000]: monthly precipitation (mm)
  stats.json         global, places and canon checks
  preview.png        annual precipitation, equirectangular
Deterministic.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
sys.path.insert(0, str(WORLD / "climate" / "scripts"))
sys.path.insert(0, str(HERE))
import cubesphere as cs     # noqa: E402
import ebm                  # noqa: E402
import climate as tlayer    # noqa: E402  (land fraction and cube interpolation helpers of the temperature layer)
import moisture             # noqa: E402

PLACES = tlayer.PLACES
DESERTS = {"11": "Dornish desert", "23": "Great Sand Sea", "24": "Red waste", "27": "Grey Waste and Cannibal Sands"}


def coarse_mean(field_cube, res, sub=4):
    lat, lon = moisture.grid(res)
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros((len(lat), len(lon)))
    for a in off:
        for b in off:
            LA, LO = np.meshgrid(np.clip(lat + a * res, -89.99, 89.99), lon + b * res, indexing="ij")
            f, j, i = cs.latlon_to_cell(np.radians(LA), np.radians(LO), cs.N)
            acc += field_cube[f, j, i]
    return acc / sub ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--relief", required=True)
    ap.add_argument("--canon-relief", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "precipitation.json").read_text())
    tcfg = json.loads((WORLD / "climate" / "data" / "climate.json").read_text())
    elev = np.load(args.relief)["elevation_m"]
    fid = np.load(args.canon_relief)["feature_id"]
    land = elev > 0

    # temperature: the same model run as the temperature layer (sea level, 2 degrees), then to 1 degree
    p = tcfg["ebm"]
    p["qflux"] = dict(p["qflux"], lon_min=tcfg["world_qflux_box"]["lon_min"], lon_max=tcfg["world_qflux_box"]["lon_max"])
    lf2 = tlayer.land_fraction(land, p["res_deg"])
    lat2, lon2, T2, ice2 = ebm.run(lf2, p, tcfg["orbit"])
    res = cfg["precip"]["res_deg"]
    lat, lon = moisture.grid(res)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    T = np.stack([tlayer.to_cube(T2[m], lat2, lon2, LA, LO) for m in range(12)])
    ice = np.stack([tlayer.to_cube(ice2[m].astype(float), lat2, lon2, LA, LO) > 0.5 for m in range(12)])
    lf = coarse_mean(land.astype(float), res)
    z = coarse_mean(np.where(land, elev, 0.0), res)
    P, E, _, diag = moisture.run(lat, lon, T, lf, z, ice, cfg)

    n = cs.N
    clat = np.zeros((6, n, n))
    clon = np.zeros((6, n, n))
    for f in range(6):
        la, lo = cs.cell_latlon(f, n)
        clat[f], clon[f] = np.degrees(la), np.degrees(lo)
    pm = np.zeros((12, 6, n, n), np.uint16)
    for m in range(12):
        pm[m] = np.clip(np.round(tlayer.to_cube(P[m], lat, lon, clat, clon)), 0, 65535).astype(np.uint16)
    np.savez_compressed(out / "precipitation.npz", p_month_mm=pm)

    A = np.broadcast_to(cs.cell_area_km2(n), clat.shape)
    ann = pm.astype(np.float32).sum(0)

    def at(la, lo):
        f, j, i = cs.latlon_to_cell(np.radians(la), np.radians(lo), n)
        return {"year_mm": int(ann[f, j, i]), "jan_mm": int(pm[0, f, j, i]), "jul_mm": int(pm[6, f, j, i])}

    deserts = {}
    for k, name in DESERTS.items():
        m = land & (fid == int(k))
        if m.any():
            deserts[name] = int(round(float((ann[m] * A[m]).sum() / A[m].sum())))
    westeros = land & (clat > 28) & (clat < 64) & (clon > -30) & (clon < 15)
    stats = {
        "global_mean_mm": int(round(float((ann * A).sum() / A.sum()))),
        "land_mean_mm": int(round(float((ann[land] * A[land]).sum() / A[land].sum()))),
        "places": {k: at(*v) for k, v in PLACES.items()},
        "desert_means_mm": deserts,
        "westeros_share_below_250mm": round(float(A[westeros & (ann < 250)].sum() / A[westeros].sum()), 3),
        "earth_reference": {"global_mean_mm": 990, "desert_mm_below": 250},
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, H = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(H) + 0.5) * 180 / H
    PLO, PLA = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(PLA), np.radians(PLO), n)
    v = np.log10(np.maximum(ann[ff, jj, ii], 10))
    stops = np.array([[1.0, 150, 60, 20], [2.0, 230, 190, 90], [2.5, 240, 240, 150], [2.8, 150, 210, 110],
                      [3.1, 40, 150, 70], [3.4, 20, 110, 160], [3.8, 40, 40, 140]], float)
    img = np.stack([np.interp(v, stops[:, 0], stops[:, c + 1]) for c in range(3)], -1)
    lnd = land[ff, jj, ii]
    img[~lnd] = img[~lnd] * 0.55 + 255 * 0.45 * np.array([0.8, 0.9, 1.0])
    edge = (lnd ^ np.roll(lnd, 1, 0)) | (lnd ^ np.roll(lnd, 1, 1))
    img[edge] = 0
    Image.fromarray(img.astype(np.uint8)).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
