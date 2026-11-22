"""Layer "vegetation": the potential natural vegetation of every cell (experiment E15).

The biome of each cell from the bioclimatic key (vegkey.py; BIOME family of models, thresholds tuned on the Earth
test): winter cold, summer heat and the moisture index (actual / potential evapotranspiration of a 150 mm bucket)
admit plant types; marshes and peat soils give wetlands (swamp forest where the climate carries forest), large rivers
through dry land give riparian woodland, salt flats stay bare. Net primary productivity from the Miami model.
"Potential" means without people: the game's clearing and farming start from this.

Input:
  --relief, --ice, --temperature, --precipitation, --hydrology, --soils (npz files of those layers)
  world/vegetation/data/vegetation.json
Output (in --out):
  vegetation.npz [6, 1000, 1000]: biome uint8 (index into "biomes"), npp_gm2 uint16 (g dry matter / m2 / yr)
  stats.json     areas, shares by region, Earth-equivalent (WWF) shares, canon checks
  preview.png    biomes, equirectangular
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
sys.path.insert(0, str(WORLD / "canon_outline" / "scripts"))
sys.path.insert(0, str(HERE))
import cubesphere as cs     # noqa: E402
import flow                 # noqa: E402
import vegkey               # noqa: E402
import world_climate as wc  # noqa: E402

REGIONS = {"Westeros south of the Wall": (28.0, 63.75, -30.0, 15.0), "Beyond the Wall": (63.75, 90.0, -40.0, 30.0),
           "Western Essos": (25.0, 55.0, 15.0, 50.0)}
# canon checks: the Haunted Forest box of the climate layer (E09), the Dothraki sea's box on the canon map (E02,
# projection study), and the canon features of world/canon_relief
DOTHRAKI_SEA_PX = (3830, 2680, 6120, 3830)


def main():
    ap = argparse.ArgumentParser()
    for a in ("--relief", "--ice", "--temperature", "--precipitation", "--hydrology", "--soils", "--canon-relief",
              "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "vegetation.json").read_text())
    scfg = json.loads((WORLD / "soils" / "data" / "soils.json").read_text())
    n = cs.N
    elev = np.load(args.relief)["elevation_m"].astype(np.float64)
    land = elev > 0
    ice = np.load(args.ice)
    surface = np.where(land, ice["surface_m"], 0.0)
    sheet = ice["sheet_frac"] == 255
    T = np.load(args.temperature)["t_month_dC"].astype(np.float32) / 10
    P = np.load(args.precipitation)["p_month_mm"].astype(np.float32)
    hy = np.load(args.hydrology)
    rsg = np.load(args.soils)["rsg"]
    clat, clon = wc.cell_latlon_deg(n)
    pet = flow.pet_thornthwaite(T, clat.astype(np.float32))
    inp = {"T": T, "P": P, "PET": pet, "elev_m": surface, "lat": clat, "land": land, "ice": sheet,
           "lake": hy["lake"], "q_m3s": hy["discharge_m3s"], "peat": rsg == scfg["groups"].index("Histosols")}
    code, climatic, bc = vegkey.classify(inp, cfg)
    biomes = cfg["biomes"]
    B = {b: i for i, b in enumerate(biomes)}
    veg = land & ~sheet & (code != B["Water"])
    npp = np.where(veg & (code != B["Salt flat"]), vegkey.npp_miami(bc["MAT"], bc["P_year"]), 0)
    np.savez_compressed(out / "vegetation.npz", biome=code.astype(np.uint8),
                        npp_gm2=np.round(npp).astype(np.uint16))

    A = np.broadcast_to(cs.cell_area_km2(n), land.shape)
    lnd = land & ~sheet
    area = {b: round(float(A[lnd & (code == k)].sum()) / 1e6, 2) for k, b in enumerate(biomes)}
    # Earth-equivalent shares of the vegetated land, against Earth's WWF shares
    wet = code == B["Wetland"]
    cold_wet = (wet & (bc["GDD5"] < cfg["wetland_wwf"]["boreal_gdd5"])) | (code == B["Swamp forest"])
    wwf = {}
    for b, w in cfg["wwf_map"].items():
        if w is not None:
            m = veg & ((code == B[b]) | (cold_wet & (climatic == B[b])))
            wwf[w] = wwf.get(w, 0.0) + float(A[m].sum())
    wwf["Flooded grasslands"] = float(A[veg & wet & ~cold_wet].sum())
    tot = sum(wwf.values())
    ref = {k: v for k, v in cfg["wwf_area_Mkm2"].items() if not k.startswith("_")}
    rtot = sum(ref.values())
    wwf_share = {k: (round(wwf.get(k, 0) / tot, 3), round(v / rtot, 3)) for k, v in ref.items()}

    regions = {}
    for name, (a0, a1, o0, o1) in REGIONS.items():
        m = lnd & (clat >= a0) & (clat < a1) & (clon >= o0) & (clon < o1)
        t = A[m].sum()
        shares = {b: round(float(A[m & (code == k)].sum() / t), 3) for k, b in enumerate(biomes)}
        regions[name] = dict(sorted(((b, v) for b, v in shares.items() if v >= 0.01), key=lambda x: -x[1]))
        regions[name]["mean_npp_gm2"] = round(float((npp * A)[m].sum() / t))

    def shares_in(m):
        t = A[m].sum()
        s = {b: round(float(A[m & (code == k)].sum() / t), 3) for k, b in enumerate(biomes)}
        return dict(sorted(((b, v) for b, v in s.items() if v >= 0.02), key=lambda x: -x[1]))

    import canon_outline as co
    x, y = co.map_px(np.radians(clat), np.radians(clon))
    x0, y0, x1, y1 = DOTHRAKI_SEA_PX
    dothraki = lnd & (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
    hf = lnd & (clat > 64) & (clat < 66) & (np.abs(clon - 0.5) < 4)
    feats = json.loads((WORLD / "canon_relief" / "data" / "canon_features.json").read_text())["features"]
    fid = np.load(args.canon_relief)["feature_id"]
    canon = {"Haunted Forest (64-66 N by the Wall)": shares_in(hf), "Dothraki sea (canon map box)": shares_in(dothraki)}
    for i, f in enumerate(feats):
        if f["type"] in ("desert", "wetland"):
            canon[f["name"]] = shares_in(lnd & (fid == i + 1))
    import climate as tlayer
    places = {}
    for name, (la, lo) in tlayer.PLACES.items():
        f, j, i = cs.latlon_to_cell(np.radians(la), np.radians(lo), n)
        places[name] = biomes[code[f, j, i]]
    stats = {"area_Mkm2": area, "wwf_share_world_vs_earth": wwf_share, "regions": regions, "canon": canon,
             "places": places, "npp_Pg_dm": round(float((npp * A * 1e6).sum()) / 1e15, 1)}
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, H = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(H) + 0.5) * 180 / H
    PLO, PLA = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(PLA), np.radians(PLO), n)
    pal = np.array([cfg["colours"][b] for b in biomes], np.uint8)
    Image.fromarray(pal[code[ff, jj, ii]]).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
