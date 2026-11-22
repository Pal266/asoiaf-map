"""Layer "soils": the dominant WRB Reference Soil Group and a natural fertility index of every cell (experiment E13).

Built from the climate with the ice, the relief, the ice and the hydrology, by the key in soilkey.py:
- monthly temperature and precipitation (temperature and precipitation layers run with --ice);
- Thornthwaite potential evapotranspiration and moisture index;
- local relief: height range over the 5 x 5 cells around each cell (~50 km, as the Earth test's 3 x 3 10' cells),
  mapped by rank onto the Earth test's distribution;
- hydrology: lakes, marshes, salt flats, closed basins and discharge;
- parent material: the lithology layer (E14: rock classes and surface deposits), and volcanic rock where the canon
  map draws volcanic relief.

Input:
  --relief, --ice, --temperature, --precipitation, --hydrology, --canon-relief, --lithology (npz files of those
  layers)
  world/soils/data/soils.json, world/canon_relief/data/canon_features.json
Output (in --out):
  soils.npz   [6, 1000, 1000]: rsg uint8 (index into soils.json "groups"), fertility uint8 (0-100)
  stats.json  areas per group, by region, places
  preview.png soil groups, equirectangular
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
sys.path.insert(0, str(WORLD / "hydrology" / "scripts"))
sys.path.insert(0, str(HERE))
import cubesphere as cs     # noqa: E402
import world_climate as wc  # noqa: E402
import flow                 # noqa: E402
import soilkey              # noqa: E402

COLOURS = {"Ice": (245, 248, 252), "Water": (60, 110, 190), "Histosols": (90, 60, 50), "Solonchaks": (230, 225, 235),
           "Cryosols": (170, 190, 210), "Leptosols": (160, 160, 160), "Andosols": (200, 80, 60),
           "Fluvisols": (110, 170, 200), "Gleysols": (100, 140, 150), "Arenosols": (245, 220, 150),
           "Gypsisols": (250, 240, 200), "Calcisols": (240, 200, 140), "Kastanozems": (190, 150, 100),
           "Chernozems": (60, 50, 40), "Phaeozems": (120, 100, 70), "Vertisols": (110, 90, 110), "Nitisols": (170, 60, 50),
           "Ferralsols": (200, 70, 40), "Acrisols": (230, 150, 90), "Lixisols": (230, 180, 120),
           "Luvisols": (200, 130, 110), "Albeluvisols": (210, 190, 170), "Podzols": (170, 160, 190),
           "Cambisols": (220, 200, 120), "Umbrisols": (120, 110, 90)}
REGIONS = {"Westeros south of the Wall": (28.0, 63.75, -30.0, 15.0), "Beyond the Wall": (63.75, 90.0, -40.0, 30.0),
           "Western Essos": (25.0, 55.0, 15.0, 50.0)}


def local_relief(elev, land, n, half=2):
    F, J, I = cs.halo_index(n, half)
    z = np.where(land, elev, 0.0)
    pad = z[F, J, I]
    k = 2 * half + 1
    mx = np.full(z.shape, -1e9)
    mn = np.full(z.shape, 1e9)
    for a in range(k):
        for b in range(k):
            win = pad[:, a:a + n, b:b + n]
            mx = np.maximum(mx, win)
            mn = np.minimum(mn, win)
    return mx - mn


def main():
    ap = argparse.ArgumentParser()
    for a in ("--relief", "--ice", "--temperature", "--precipitation", "--hydrology", "--canon-relief", "--out"):
        ap.add_argument(a, required=True)
    ap.add_argument("--lithology", default=None, help="lithology.npz from world/lithology/scripts/lithology.py")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "soils.json").read_text())
    n = cs.N
    elev = np.load(args.relief)["elevation_m"].astype(np.float64)
    land = elev > 0
    ice = np.load(args.ice)
    surface = np.where(land, ice["surface_m"], 0.0)
    sheet = ice["sheet_frac"] == 255
    T = np.load(args.temperature)["t_month_dC"].astype(np.float32) / 10
    P = np.load(args.precipitation)["p_month_mm"].astype(np.float32)
    hy = np.load(args.hydrology)
    clat, clon = wc.cell_latlon_deg(n)
    pet = flow.pet_thornthwaite(T, clat.astype(np.float32))
    feats = json.loads((WORLD / "canon_relief" / "data" / "canon_features.json").read_text())["features"]
    volc_ids = [i + 1 for i, f in enumerate(feats) if f["type"] == "volcanic"]
    fid = np.load(args.canon_relief)["feature_id"]
    # local relief, mapped by rank onto the Earth test's distribution (the key's thresholds are set on it)
    rel = local_relief(surface, land, n)
    A = np.broadcast_to(cs.cell_area_km2(n), land.shape)
    m = land & ~sheet
    order = np.argsort(rel[m])
    cum = np.cumsum(A[m][order])
    rank = np.empty(order.size)
    rank[order] = cum / cum[-1]
    rel_q = np.zeros(rel.shape)
    rel_q[m] = np.interp(rank, np.linspace(0, 1, 101), cfg["earth_relief_quantiles_m"])
    inp = {"T": T, "P": P, "PET": pet, "relief_m": rel_q, "elev_m": surface,
           "land": land, "ice": sheet, "lake": hy["lake"], "q_m3s": hy["discharge_m3s"],
           "volcanic": land & np.isin(fid, volc_ids), "basin": np.isin(hy["lake"], (2, 3, 4))}
    if args.lithology:
        li = np.load(args.lithology)
        inp.update(rock=li["rock"], deposit=li["deposit"], bedrock=li["bedrock"])
    code = soilkey.classify(inp, cfg)
    groups = cfg["groups"]
    fert = np.array([cfg["fertility"][g] for g in groups])[code]
    np.savez_compressed(out / "soils.npz", rsg=code.astype(np.uint8), fertility=np.round(fert * 100).astype(np.uint8))

    lnd = land & ~sheet
    area = {g: round(float(A[lnd & (code == k)].sum()) / 1e4) for k, g in enumerate(groups)}
    regions = {}
    for name, (a0, a1, o0, o1) in REGIONS.items():
        m = lnd & (clat >= a0) & (clat < a1) & (clon >= o0) & (clon < o1)
        tot = A[m].sum()
        shares = {g: round(float(A[m & (code == k)].sum() / tot), 3) for k, g in enumerate(groups)}
        regions[name] = dict(sorted(((g, v) for g, v in shares.items() if v >= 0.01), key=lambda x: -x[1]))
        regions[name]["mean_fertility"] = round(float((fert * A)[m].sum() / tot), 2)
    import climate as tlayer
    places = {}
    for name, (la, lo) in tlayer.PLACES.items():
        f, j, i = cs.latlon_to_cell(np.radians(la), np.radians(lo), n)
        places[name] = groups[code[f, j, i]]
    stats = {"area_Mha": area, "regions": regions, "places": places, "fao_earth_area_Mha": cfg["fao_area_Mha"]}
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, H = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(H) + 0.5) * 180 / H
    PLO, PLA = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(PLA), np.radians(PLO), n)
    pal = np.array([COLOURS[g] for g in groups], np.uint8)
    Image.fromarray(pal[code[ff, jj, ii]]).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
