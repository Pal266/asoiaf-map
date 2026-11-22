"""Test the soil key on Earth (experiment E13).

Earth inputs as for the ice and hydrology tests: temperature and precipitation of the E09 and E10 models at 10',
ETOPO1 surface and today's ice (world/ice/scripts/earth_ice_validation.py), and the Earth hydrology (lakes, marshes,
salt flats, discharge) of world/hydrology/scripts/earth_hydro_validation.py. Local relief: height range in the
3 x 3 neighbourhood (~55 km). No volcanic parent material is given, so Andosols are not tested.

Scores:
- sites: the share of data/earth_soil_sites.json whose group is among the accepted ones (at the site's cell or
  one of its 8 neighbours);
- areas: modelled area of each group against the FAO estimates (data/soils.json, fao_area_Mha), log10 RMSE.

Usage: python3 earth_soil_validation.py --earth-data DIR
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
sys.path.insert(0, str(WORLD / "hydrology" / "scripts"))
sys.path.insert(0, str(HERE))
import earth_ice_validation as iv     # noqa: E402
import earth_hydro_validation as hv   # noqa: E402
import flow                           # noqa: E402
import soilkey                        # noqa: E402

RES = 1 / 6


def earth_inputs(earth, hcfg, nbd=None):
    """Soil-key inputs on Earth's 10' grid from the ice-test climate and the hydrology test."""
    lat, lon = earth["lat"], earth["lon"]
    ny, nx = len(lat), len(lon)
    _, res = hv.evaluate(hcfg, earth, nbd)
    lake = np.select([res["terminal"], res["lake"], res["marsh"], res["playa"]], [2, 1, 3, 4], 0).reshape(ny, nx)
    LAT = np.repeat(lat, nx).reshape(ny, nx)
    pet = flow.pet_thornthwaite(earth["T"].astype(np.float64), LAT)      # plain Thornthwaite, as his classes
    z = np.where(earth["land"], earth["z"], 0.0)
    zz = np.pad(z, ((1, 1), (0, 0)), mode="edge")
    stack = [np.roll(zz, d, 1)[1 + a:1 + a + ny] for a in (-1, 0, 1) for d in (-1, 0, 1)]
    relief = np.max(stack, 0) - np.min(stack, 0)
    return {"T": earth["T"], "P": earth["P"], "PET": pet, "relief_m": relief, "elev_m": z,
            "land": earth["land"], "ice": earth["obs"], "lake": lake,
            "q_m3s": res["discharge_m3s"].reshape(ny, nx), "volcanic": np.zeros((ny, nx), bool),
            "basin": (res["filled"] - (earth["z"] - hcfg["earth_test"]["burn_m"] * earth["rivers"]).ravel()
                      > 1).reshape(ny, nx)}


def evaluate(cfg, inp, lat, lon):
    code = soilkey.classify(inp, cfg)
    groups = cfg["groups"]
    ny, nx = code.shape
    A = (np.cos(np.radians(lat)) * (111.2 * RES) ** 2)[:, None] * np.ones((1, nx))
    sites = json.loads((DATA / "earth_soil_sites.json").read_text())["sites"]
    rows, hits = [], 0
    for s in sites:
        j = int((90 - s["lat"]) / RES)
        i = int((s["lon"] + 180) / RES) % nx
        win = [groups[code[jj, ii % nx]] for jj in (j - 1, j, j + 1) for ii in (i - 1, i, i + 1) if 0 <= jj < ny]
        here = groups[code[j, i]]
        ok = any(g in s["accept"] for g in win)
        hits += ok
        rows.append([s["name"], here, ok])
    area = {g: round(float(A[code == k].sum()) / 1e4) for k, g in enumerate(groups)}      # Mha
    fao = cfg["fao_area_Mha"]
    e = [np.log10(max(area[g], 1) / fao[g]) for g in fao if g != "Andosols"]
    return {"site_hit_rate": round(hits / len(sites), 3), "area_rmse_log10": round(float(np.sqrt(np.mean(np.square(e)))), 3),
            "area_Mha": area, "sites": rows}, code


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--earth-data", required=True)
    args = ap.parse_args()
    cfg = json.loads((DATA / "soils.json").read_text())
    hcfg = json.loads((WORLD / "hydrology" / "data" / "hydrology.json").read_text())
    earth = iv.load_earth(args.earth_data)
    earth["rivers"] = hv.river_mask(args.earth_data, earth["lat"], earth["lon"]) & earth["land"]
    inp = earth_inputs(earth, hcfg)
    out, _ = evaluate(cfg, inp, earth["lat"], earth["lon"])
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
