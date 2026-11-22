"""Test the vegetation key on Earth (experiment E15).

Earth inputs as for the soil test (world/soils/scripts/earth_soil_validation.py): temperature and precipitation of
the E09 and E10 models at 10', ETOPO1 surface, today's ice, the Earth hydrology, and the soil key's Histosols on
Earth for peat. The climate is the model's, so the key's thresholds absorb the model's errors on Earth in the same
way the world run sees them.

Scores:
- sites: the share of data/earth_veg_sites.json whose biome is among the accepted ones (at the site's cell or one of
  its 8 neighbours);
- areas: modelled area of each WWF biome (Olson et al. 2001) against the tabulated areas, log10 RMSE over the
  biomes the key produces (mangroves and tropical conifer forests are not modelled);
- NPP: global Miami-model productivity against ~120 Pg dry matter per year.

Usage: python3 earth_veg_validation.py --earth-data DIR
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
DATA = HERE.parent / "data"
sys.path.insert(0, str(WORLD / "soils" / "scripts"))
sys.path.insert(0, str(HERE))
import vegkey                          # noqa: E402

RES = 1 / 6


def earth_inputs(soil_inp, soil_code, soil_groups):
    """Vegetation inputs from the soil test's Earth inputs and its soil groups."""
    lat = soil_inp["lat"]
    ny, nx = soil_inp["land"].shape
    return {"T": soil_inp["T"], "P": soil_inp["P"], "PET": soil_inp["PET"], "elev_m": soil_inp["elev_m"],
            "lat": np.repeat(lat, nx).reshape(ny, nx), "land": soil_inp["land"], "ice": soil_inp["ice"],
            "lake": soil_inp["lake"], "q_m3s": soil_inp["q_m3s"],
            "peat": soil_code == soil_groups.index("Histosols")}


def evaluate(cfg, inp, lat, lon):
    code, climatic, bc = vegkey.classify(inp, cfg)
    biomes = cfg["biomes"]
    B = {b: i for i, b in enumerate(biomes)}
    ny, nx = code.shape
    A = (np.cos(np.radians(lat)) * (111.2 * RES) ** 2)[:, None] * np.ones((1, nx))
    # WWF biome of every cell
    wwf = {}
    land = inp["land"] & ~inp["ice"] & (code != B["Water"])
    # wetlands: WWF puts peatlands of the cold zones in boreal and tundra and swamp forests in the forest biomes;
    # only open wetlands of warmer climates are "flooded grasslands"
    wet = code == B["Wetland"]
    cold_wet = (wet & (bc["GDD5"] < cfg["wetland_wwf"]["boreal_gdd5"])) | (code == B["Swamp forest"])
    for b, w in cfg["wwf_map"].items():
        if w is None:
            continue
        m = land & ((code == B[b]) | (cold_wet & (climatic == B[b])))
        wwf[w] = wwf.get(w, 0.0) + float(A[m].sum()) / 1e6
    wwf["Flooded grasslands"] = float(A[land & wet & ~cold_wet].sum()) / 1e6
    ref = {k: v for k, v in cfg["wwf_area_Mkm2"].items() if not k.startswith("_")}
    used = [k for k in ref if wwf.get(k, 0) > 0]
    rmse = float(np.sqrt(np.mean([np.log10(wwf[k] / ref[k]) ** 2 for k in used])))
    sites = json.loads((DATA / "earth_veg_sites.json").read_text())["sites"]
    res = []
    for name, la, lo, ok in sites:
        j = int(np.argmin(np.abs(lat - la)))
        i = int(np.argmin(np.abs(((lon - lo + 180) % 360) - 180)))
        near = {biomes[code[jj, ii % nx]] for jj in range(max(j - 1, 0), min(j + 2, ny)) for ii in (i - 1, i, i + 1)}
        res.append((name, biomes[code[j, i]], bool(near & set(ok)), ok))
    npp = vegkey.npp_miami(bc["MAT"], bc["P_year"])
    npp_pg = float((npp * A * 1e6)[land].sum()) / 1e15
    return {"site_hit_rate": round(float(np.mean([r[2] for r in res])), 3), "area_rmse_log10": round(rmse, 3),
            "wwf_area_Mkm2": {k: (round(wwf.get(k, 0), 2), ref[k]) for k in ref},
            "npp_Pg_dm": round(npp_pg, 1), "sites": res}, code


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--earth-data", required=True, help="directory with the Earth data of the ice and hydrology tests")
    args = ap.parse_args()
    sys.path.insert(0, str(WORLD / "ice" / "scripts"))
    sys.path.insert(0, str(WORLD / "hydrology" / "scripts"))
    import earth_hydro_validation as hv
    import earth_ice_validation as iv
    import earth_soil_validation as sv
    import soilkey
    hcfg = json.loads((WORLD / "hydrology" / "data" / "hydrology.json").read_text())
    scfg = json.loads((WORLD / "soils" / "data" / "soils.json").read_text())
    earth = iv.load_earth(args.earth_data)
    earth["rivers"] = hv.river_mask(args.earth_data, earth["lat"], earth["lon"]) & earth["land"]
    sinp = sv.earth_inputs(earth, hcfg)
    sinp["lat"], sinp["lon"] = earth["lat"], earth["lon"]
    scode = soilkey.classify(sinp, scfg)
    cfg = json.loads((DATA / "vegetation.json").read_text())
    out, _ = evaluate(cfg, earth_inputs(sinp, scode, scfg["groups"]), earth["lat"], earth["lon"])
    print(json.dumps({k: v for k, v in out.items() if k != "sites"}, indent=1))
    for r in out["sites"]:
        print(("  ok  " if r[2] else "  MISS"), r[0], r[1], r[3])


if __name__ == "__main__":
    main()
