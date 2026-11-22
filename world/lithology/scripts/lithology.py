"""Layer "lithology": the rock at the surface of every cell (experiment E14).

Rock classes are those of the Global Lithological Map (GLiM, Hartmann & Moosdorf 2012, level 1).
1. Bedrock from the tectonic province of the geological history (E07), read at a point displaced by a smooth random
   warp (~0.8°) so that geological contacts are not the history's straight polygon edges. Each province has an
   assumed composition (data/lithology.json); cratons are split between exposed shield and platform cover by local relief (shields
   crop out where erosion has stripped the cover). Within a province the class of a cell comes from one smooth
   random field (two octaves, ~130 km and ~30 km), so that rocks form coherent outcrops with the province's shares.
2. Canon volcanic relief: intermediate volcanics and pyroclastics.
3. Surface deposits over the bedrock (GLiM "unconsolidated sediments"), from the hydrology and climate:
   alluvium along large rivers and in marshes, evaporites on salt flats, aeolian sand seas on low desert plains,
   loess on low semi-arid to subhumid plains within reach of the deserts' dust.
4. Lakes are water, ice sheets ice.

Input:
  --history, --relief, --ice, --hydrology, --precipitation, --temperature, --canon-relief (npz of those layers)
  world/lithology/data/lithology.json
Output (in --out):
  lithology.npz  [6, 1000, 1000]: rock uint8 (index into "order"), deposit uint8 (0 none, 1 alluvium,
                 2 evaporite, 3 aeolian sand, 4 loess), bedrock uint8 (the rock under the deposits)
  stats.json     land shares against GLiM, by region
  preview.png    rock classes, equirectangular
Deterministic.
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
sys.path.insert(0, str(WORLD / "climate" / "scripts"))
sys.path.insert(0, str(WORLD / "soils" / "scripts"))
import cubesphere as cs     # noqa: E402
import world_climate as wc  # noqa: E402
from soils import local_relief  # noqa: E402

COLOURS = {"su": (240, 230, 180), "ss": (230, 190, 120), "sm": (200, 170, 130), "sc": (150, 200, 230),
           "ev": (250, 200, 230), "mt": (170, 120, 170), "pa": (240, 130, 130), "pi": (220, 100, 90),
           "pb": (120, 60, 60), "va": (250, 170, 90), "vi": (230, 120, 40), "vb": (90, 90, 90),
           "py": (200, 150, 60), "ig": (245, 248, 252), "wb": (60, 110, 190)}
DEPOSITS = ["none", "alluvium", "evaporite", "aeolian sand", "loess"]
REGIONS = {"Westeros south of the Wall": (28.0, 63.75, -30.0, 15.0), "Beyond the Wall": (63.75, 90.0, -40.0, 30.0),
           "Western Essos": (25.0, 55.0, 15.0, 50.0)}


def smooth_field(n_coarse, iters, n, rng, clat_rad=None):
    """A smooth standard-normal random field on the cube of resolution n, from white noise on a coarse cube
    (seam-aware box smoothing), bilinear within faces."""
    f = cs.box_smooth(rng.standard_normal((6, n_coarse, n_coarse)), iters)
    f = (f - f.mean()) / f.std()
    k = n / n_coarse
    x = (np.arange(n) + 0.5) / k - 0.5
    x0 = np.clip(np.floor(x).astype(int), 0, n_coarse - 2)
    t = np.clip(x - x0, 0, 1)
    J0, I0 = np.meshgrid(x0, x0, indexing="ij")
    TJ, TI = np.meshgrid(t, t, indexing="ij")
    out = np.zeros((6, n, n))
    for fc in range(6):
        g = f[fc]
        out[fc] = ((1 - TJ) * (1 - TI) * g[J0, I0] + (1 - TJ) * TI * g[J0, I0 + 1]
                   + TJ * (1 - TI) * g[J0 + 1, I0] + TJ * TI * g[J0 + 1, I0 + 1])
    return out


def uniform_field(cfg, n, rng):
    nz = cfg["noise"]
    a = smooth_field(nz["coarse_n"], nz["coarse_smooth"], n, rng)
    b = smooth_field(nz["fine_n"], nz["fine_smooth"], n, rng)
    g = (a + nz["fine_weight"] * b) / math.sqrt(1 + nz["fine_weight"] ** 2)
    from scipy.special import ndtr
    return ndtr(g)          # uniform on [0, 1], spatially coherent


def pick(U, comp, order):
    """Class index by the cumulative shares of comp (dict class -> share) at uniform values U."""
    names = [c for c in order if c in comp]
    p = np.array([comp[c] for c in names], float)
    cum = np.cumsum(p / p.sum())
    idx = np.minimum(np.searchsorted(cum, U), len(names) - 1)
    return np.array([order.index(c) for c in names])[idx]


def main():
    ap = argparse.ArgumentParser()
    for a in ("--history", "--relief", "--ice", "--hydrology", "--precipitation", "--temperature", "--canon-relief",
              "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "lithology.json").read_text())
    order = cfg["order"]
    O = {c: i for i, c in enumerate(order)}
    n = cs.N
    rng = np.random.default_rng(cfg["noise"]["seed"])
    prov = np.load(args.history)["province"]
    elev = np.load(args.relief)["elevation_m"].astype(np.float64)
    land = elev > 0
    sheet = np.load(args.ice)["sheet_frac"] == 255
    hy = np.load(args.hydrology)
    P = np.load(args.precipitation)["p_month_mm"].astype(np.float32).sum(0)
    T = np.load(args.temperature)["t_month_dC"].astype(np.float32) / 10
    rel = local_relief(elev, land, n)
    U = uniform_field(cfg, n, rng)
    U2 = uniform_field(cfg, n, rng)
    clat, clon = wc.cell_latlon_deg(n)
    # the provinces of the history are drawn as polygons and buffers; their boundaries are warped by a smooth random
    # displacement (geological contacts are not straight), and the province is read at the displaced point
    nz = cfg["noise"]
    amp = cfg["province_warp_deg"]
    wy = smooth_field(nz["coarse_n"], nz["coarse_smooth"], n, rng) * amp
    wx = smooth_field(nz["coarse_n"], nz["coarse_smooth"], n, rng) * amp / np.maximum(np.cos(np.radians(clat)), 0.2)
    fw, jw, iw = cs.latlon_to_cell(np.radians(np.clip(clat + wy, -89.9, 89.9)), np.radians(clon + wx), n)
    pw = prov[fw, jw, iw]
    prov = np.where(land & (pw > 0), pw, prov)

    # 1. bedrock by province; cratons: shield where relief is high, platform cover where it is low
    bed = np.full(land.shape, O["ss"], np.int16)
    comp = cfg["province_composition"]
    for code in range(1, 10):
        m = land & (prov == code)
        if not m.any():
            continue
        if code == 1:
            w = np.clip(rel[m] / (2 * cfg["craton_shield_relief_m"]), cfg["craton_shield_min_share"], 1.0)
            shield = U2[m] < w          # a second field decides shield or cover, so both keep coherent outcrops
            b = pick(U[m], comp["2"], order)
            b[shield] = pick(U[m][shield], comp["shield"], order)
            bed[m] = b
        else:
            bed[m] = pick(U[m], comp[str(code)], order)
    # 2. canon volcanic relief
    feats = json.loads((WORLD / "canon_relief" / "data" / "canon_features.json").read_text())["features"]
    volc = land & np.isin(np.load(args.canon_relief)["feature_id"],
                          [i + 1 for i, f in enumerate(feats) if f["type"] == "volcanic"])
    bed[volc] = np.where(U[volc] < 0.75, O["vi"], O["py"])

    # 3. surface deposits
    dp = cfg["deposits"]
    lake = hy["lake"]
    dep = np.zeros(land.shape, np.uint8)
    alluv = land & (((hy["discharge_m3s"] >= dp["alluvium_q_m3s"]) & (rel < dp["alluvium_relief_m"])) | (lake == 3))
    dep[alluv] = 1
    dep[land & (lake == 4)] = 2
    # warm deserts only: polar deserts are dry but make frost debris, not sand seas or dust
    arid = land & (P < dp["aeolian_p_mm"]) & (T.mean(0) >= dp["aeolian_mat_min_C"])
    dep[(dep == 0) & arid & (rel < dp["aeolian_relief_m"]) & (U2 < dp["erg_share"])] = 3
    # loess: dust from the deserts settles on the semi-arid and subhumid plains around them (~800 km)
    nc = cfg["noise"]["coarse_n"]
    fc, jc, ic = cs.latlon_to_cell(np.radians(clat), np.radians(clon), nc)
    src = np.zeros((6, nc, nc))
    np.add.at(src, (fc[arid], jc[arid], ic[arid]), 1.0)
    cnt = np.zeros((6, nc, nc))
    np.add.at(cnt, (fc, jc, ic), 1.0)
    src = src / np.maximum(cnt, 1)
    cell_km = 2 * math.pi * cs.R_KM / (4 * nc)
    sigma = dp["loess_distance_km"] / cell_km / 2
    near = cs.box_smooth(src, max(1, int(round(1.5 * sigma ** 2))))[fc, jc, ic]
    loess = (land & (dep == 0) & ~arid & (P < 800) & (near > 0.02) & (rel < dp["loess_relief_m"]) &
             (T.min(0) < dp["loess_cold_month_max_C"]))
    dep[loess] = 4

    rock = bed.copy()
    rock[dep > 0] = np.where(dep[dep > 0] == 2, O["ev"], O["su"])
    rock[sheet] = O["ig"]
    rock[land & np.isin(lake, (1, 2))] = O["wb"]
    rock[~land] = O["wb"]
    np.savez_compressed(out / "lithology.npz", rock=rock.astype(np.uint8), deposit=dep, bedrock=bed.astype(np.uint8))

    A = np.broadcast_to(cs.cell_area_km2(n), land.shape)
    L = land
    tot = float(A[L].sum())
    share = {c: round(float(A[L & (rock == k)].sum()) / tot, 3) for k, c in enumerate(order)}
    groups = {"sediments": ["su", "ss", "sm", "sc", "ev"], "metamorphic": ["mt"], "plutonic": ["pa", "pi", "pb"],
              "volcanic": ["va", "vi", "vb", "py"], "water_ice": ["ig", "wb"]}
    gshare = {g: round(sum(share[c] for c in cs_), 3) for g, cs_ in groups.items()}
    gshare["carbonate_share_of_sediments"] = round(share["sc"] / max(gshare["sediments"], 1e-9), 3)
    gshare["unconsolidated"] = share["su"]
    deps = {d: round(float(A[L & (dep == k)].sum()) / tot, 3) for k, d in enumerate(DEPOSITS) if k}
    regions = {}
    for name, (a0, a1, o0, o1) in REGIONS.items():
        m = L & (clat >= a0) & (clat < a1) & (clon >= o0) & (clon < o1)
        t = A[m].sum()
        regions[name] = {c: round(float(A[m & (rock == k)].sum() / t), 3) for k, c in enumerate(order)
                         if A[m & (rock == k)].sum() / t >= 0.01}
    stats = {"land_shares": share, "groups": gshare, "glim": cfg["glim_land_shares"], "deposits": deps,
             "regions": regions}
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, H = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(H) + 0.5) * 180 / H
    PLO, PLA = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(PLA), np.radians(PLO), n)
    pal = np.array([COLOURS[c] for c in order], np.uint8)
    Image.fromarray(pal[rock[ff, jj, ii]]).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
