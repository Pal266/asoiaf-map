"""Layer "climate: temperature": monthly mean near-surface air temperature of every cell of the global cubed sphere
(experiment E09).

The seasonal energy balance model (ebm.py, calibrated on Earth by earth_validation.py) runs on the world's land
mask at the coarse resolution of climate.json, giving sea-level temperatures; these are interpolated to the 10 km
grid and reduced to the ground with the lapse rate.

Input:
  --relief   relief.npz (elevation_m) from world/relief/scripts/relief.py
  world/climate/data/climate.json (in the repository)
Output (in --out):
  temperature.npz  t_month_dC int16 [12, 6, 1000, 1000]: monthly mean temperature in tenths of a degree C
                   (ground level on land, sea surface at sea); sea_ice uint8 [12, ny, nx] on the coarse grid
  stats.json       global and canon checks
  preview.png      warmest and coldest month, equirectangular
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
import cubesphere as cs     # noqa: E402
import ebm                  # noqa: E402
import world_climate as wc  # noqa: E402

CFG = HERE.parent / "data" / "climate.json"
PLACES = {"Castle Black": (63.77, 1.26), "Winterfell": (56.51, -4.25), "King's Landing": (38.69, 1.18),
          "Oldtown": (29.99, -11.09), "Sunspear": (28.74, 4.99), "Volantis": (28.74, 22.86),
          "Qarth": (22.46, 57.32), "Asshai": (9.03, 83.02)}


def land_fraction(land_cube, res, sub=4):
    lat, lon = ebm.grid(res)
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros((len(lat), len(lon)))
    for a in off:
        for b in off:
            LA, LO = np.meshgrid(np.clip(lat + a * res, -89.99, 89.99), lon + b * res, indexing="ij")
            f, j, i = cs.latlon_to_cell(np.radians(LA), np.radians(LO), cs.N)
            acc += land_cube[f, j, i]
    return acc / sub ** 2


def to_cube(field, lat, lon, cell_lat, cell_lon):
    """Bilinear interpolation of a coarse [ny, nx] field (lat descending, periodic lon) to cell centres."""
    res = lat[0] - lat[1]
    y = np.clip((lat[0] - cell_lat) / res, 0, len(lat) - 1.000001)
    x = (cell_lon - lon[0]) / res
    j0 = np.floor(y).astype(int)
    fy = y - j0
    j1 = np.minimum(j0 + 1, len(lat) - 1)
    i0 = np.floor(x).astype(int) % len(lon)
    fx = x - np.floor(x)
    i1 = (i0 + 1) % len(lon)
    return ((1 - fy) * ((1 - fx) * field[j0, i0] + fx * field[j0, i1])
            + fy * ((1 - fx) * field[j1, i0] + fx * field[j1, i1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--relief", required=True)
    ap.add_argument("--ice", default=None, help="ice.npz from world/ice/scripts/ice.py (optional)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    elev = np.load(args.relief)["elevation_m"]
    land = elev > 0
    ground = np.where(land, elev, 0.0)
    sheet = None
    if args.ice:
        # with the ice layer: the ground is the ice surface, and ice sheets keep the ice albedo all year
        ice = np.load(args.ice)
        ground = np.where(land, ice["surface_m"], 0.0)
        sheet = ice["sheet_frac"].astype(np.float32) / 255
    n = cs.N
    clat, clon = wc.cell_latlon_deg(n)
    temp = wc.temperature(ground, land, sheet, clat, clon)
    t = np.round(temp["t_cube"] * 10).astype(np.int16)
    ice = temp["sea_ice"]
    np.savez_compressed(out / "temperature.npz", t_month_dC=t, sea_ice=ice.astype(np.uint8))

    # ---- stats and canon checks ----
    A = np.broadcast_to(cs.cell_area_km2(n), clat.shape)
    tm = t.astype(np.float32) / 10
    ann = tm.mean(0)
    warm = tm.max(0)
    cold = tm.min(0)

    def at(la, lo):
        f, j, i = cs.latlon_to_cell(np.radians(la), np.radians(lo), n)
        return {"elev_m": round(float(elev[f, j, i])), "jan": round(float(tm[0, f, j, i]), 1),
                "jul": round(float(tm[6, f, j, i]), 1), "year": round(float(ann[f, j, i]), 1)}

    # the Haunted Forest: land just north of the Wall (64-66 N, the Wall's longitudes +-4 deg)
    hf = land & (clat > 64) & (clat < 66) & (np.abs(clon - 0.5) < 4)
    law = land & (clat > 75) & (clon > -45) & (clon < 5)
    stats = {
        "global_mean_C": round(float((ann * A).sum() / A.sum()), 2),
        "land_mean_C": round(float((ann[land] * A[land]).sum() / A[land].sum()), 2),
        "places": {k: at(*v) for k, v in PLACES.items()},
        "haunted_forest": {"warmest_month_mean_C": round(float((warm[hf] * A[hf]).sum() / A[hf].sum()), 1),
                           "share_with_warmest_month_ge_10C": round(float(A[hf & (warm >= 10)].sum() / A[hf].sum()), 3)},
        "land_of_always_winter": {"warmest_month_mean_C": round(float((warm[law] * A[law]).sum() / A[law].sum()), 1),
                                  "share_with_warmest_month_lt_0C": round(float(A[law & (warm < 0)].sum() / A[law].sum()), 3)},
        "earth_reference": {"global_mean_C": 14.0, "tree_line_warmest_month_C": 10.0},
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    # ---- preview: July (top) and January (bottom) ----
    W, H = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(H) + 0.5) * 180 / H
    LON, LAT = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(LAT), np.radians(LON), n)
    stops = np.array([[-40, 120, 0, 160], [-20, 40, 60, 220], [0, 220, 240, 255], [10, 120, 200, 120],
                      [20, 240, 220, 80], [30, 230, 90, 40], [40, 120, 0, 0]], float)
    panels = []
    for m in (6, 0):
        v = tm[m][ff, jj, ii]
        img = np.stack([np.interp(v, stops[:, 0], stops[:, c + 1]) for c in range(3)], -1)
        lnd = land[ff, jj, ii]
        edge = (lnd ^ np.roll(lnd, 1, 0)) | (lnd ^ np.roll(lnd, 1, 1))
        img[edge] = 0
        panels.append(img)
    Image.fromarray(np.concatenate(panels, 0).astype(np.uint8)).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
