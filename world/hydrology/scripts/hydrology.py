"""Layer "hydrology": runoff, rivers and lakes of every cell of the global cubed sphere (experiment E12).

Built from the climate with the ice (temperature and precipitation layers run with --ice) and the ice surface:
1. annual potential evapotranspiration (Thornthwaite, x pet_factor) and runoff from the Budyko curve (flow.py);
   on ice sheets all precipitation runs off;
2. routing on the ice surface (bedrock where there is no ice), with the canon rivers of the map carved in
   (falling from source to mouth): priority-flood depression filling from the sea, steepest descent over 8
   neighbours, downstream accumulation;
3. lakes: depressions that overflow keep water only where deeper than lake_overflow_depth_m; depressions whose
   inflow cannot feed their open-water evaporation are terminal lakes (shrunk to the balance) and end their rivers;
   the canon lakes of the map are lakes.

Input:
  --relief         relief.npz (elevation_m)                    world/relief/scripts/relief.py
  --ice            ice.npz (surface_m, sheet_frac)              world/ice/scripts/ice.py
  --temperature    temperature.npz (t_month_dC)                 world/climate/scripts/climate.py --ice
  --precipitation  precipitation.npz (p_month_mm)               world/precipitation/scripts/precipitation.py --ice
  --canon-relief   canon_relief.npz (feature_id)                for the canon checks
  --canon-rivers   canon_rivers.npz                             world/hydrology/scripts/canon_rivers.py
  world/hydrology/data/hydrology.json
Output (in --out):
  hydrology.npz  [6, 1000, 1000]: discharge_m3s float32, runoff_mm float32, downstream int8 (index 0-7 of the
                 receiving neighbour, -1 at sea and in terminal lakes' sinks), lake uint8 (1 lake, 2 terminal lake, 3 marsh or
                 floodplain: the sediment-filled part of an overflowing depression, 4 salt flat: the dry floor of a
                 terminal depression),
                 strahler uint8 (stream order of channels with discharge >= channel_min_m3s)
  stats.json     rivers, lakes, canon checks
  preview.png    rivers and lakes over the relief
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
import world_climate as wc  # noqa: E402
import flow                 # noqa: E402

OFFSETS = ((0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1))   # (dj, di), index = downstream code


def neighbours8(n=cs.N):
    """Flat indices (8, 6 n n) of the 8 neighbours of every cell, and the distances (m) between centres."""
    F, J, I = cs.halo_index(n, 1)
    flat = (F * n * n + J * n + I).astype(np.int64)
    nb = np.stack([flat[:, 1 + dj:1 + dj + n, 1 + di:1 + di + n] for dj, di in OFFSETS]).reshape(8, -1)
    xyz = np.zeros((6 * n * n, 3))
    for f in range(6):
        la, lo = cs.cell_latlon(f, n)
        xyz[f * n * n:(f + 1) * n * n] = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                                                   np.sin(la)], -1).reshape(-1, 3)
    dist = np.stack([np.linalg.norm(xyz - xyz[nb[k]], axis=1) * cs.R_KM * 1000 for k in range(8)]) + 1e-6
    return nb, dist


def main():
    ap = argparse.ArgumentParser()
    for a in ("--relief", "--ice", "--temperature", "--precipitation", "--canon-relief", "--canon-rivers", "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "hydrology.json").read_text())
    hp = cfg["hydro"]
    n = cs.N
    elev = np.load(args.relief)["elevation_m"].astype(np.float64)
    ice = np.load(args.ice)
    land = elev > 0
    surface = np.where(land, ice["surface_m"].astype(np.float64), elev)
    sheet = ice["sheet_frac"] == 255
    T = np.load(args.temperature)["t_month_dC"].astype(np.float32) / 10
    P = np.load(args.precipitation)["p_month_mm"].astype(np.float32).sum(0)
    clat, clon = wc.cell_latlon_deg(n)
    pet = hp["pet_factor"] * flow.pet_thornthwaite(T, clat.astype(np.float32)).sum(0)
    A = np.broadcast_to(cs.cell_area_km2(n), land.shape).astype(np.float64)
    nb, dist = neighbours8(n)
    # the rivers of the official map are carved into the routing surface (decision 9): the relief (E08) has no
    # river incision, so its drainage does not follow them; lakes are found on the uncarved surface
    cr = np.load(args.canon_rivers)
    rid = cr["river_id"]
    carved = carve_canon(surface, land, cr, hp, nb)
    # the inland lakes of the map are lakes (decision 10); their water balance decides whether they overflow
    lid = cr["lake_id"].ravel()
    res = flow.route(carved.ravel(), land.ravel(), P.ravel().astype(np.float64), pet.ravel().astype(np.float64),
                     A.ravel(), nb, dist, hp, ice=sheet.ravel(), z_dep=surface.ravel(), forced_lake=lid > 0)
    lake_names = [x["name"] for x in json.loads((HERE.parent / "data" / "canon_rivers.json").read_text())["lakes"]]
    canon_lakes = {}
    for k, nm in enumerate(lake_names):
        cells = np.flatnonzero((lid == k + 1) & land.ravel())
        if cells.size == 0:
            continue
        inside = np.zeros(lid.size, bool)
        inside[cells] = True
        r_ = res["rec"][cells]
        out_cells = cells[(r_ >= 0) & ~inside[np.maximum(r_, 0)]]
        outflow = float(res["discharge_m3s"][out_cells].max()) if out_cells.size else 0.0
        if outflow < 1:
            res["terminal"][cells] = True
        canon_lakes[nm] = {"area_km2": round(float(A.ravel()[cells].sum())), "outflow_m3s": round(outflow, 1),
                           "terminal": bool(outflow < 1)}
    q = res["discharge_m3s"]
    rec = res["rec"]
    down = np.full(rec.size, -1, np.int8)
    for k in range(8):
        down[(rec >= 0) & (nb[k] == rec)] = k
    sinks = res["terminal"] & (q <= 0)
    down[sinks] = -1
    lake = np.select([res["terminal"], res["lake"], res["marsh"], res["playa"]], [2, 1, 3, 4], 0).astype(np.uint8)
    so = flow.stream_order(rec, res["order"], q, hp["channel_min_m3s"])
    shape = land.shape
    np.savez_compressed(out / "hydrology.npz", discharge_m3s=q.reshape(shape).astype(np.float32),
                        runoff_mm=res["runoff_mm"].reshape(shape).astype(np.float32),
                        downstream=down.reshape(shape), lake=lake.reshape(shape),
                        strahler=np.minimum(so, 255).astype(np.uint8).reshape(shape))

    # ---------------- stats ----------------
    landf = land.ravel()
    mouth = landf & (rec >= 0) & ~landf[np.maximum(rec, 0)]
    order_mouths = np.flatnonzero(mouth)[np.argsort(-q[mouth])]
    LA, LO = clat.ravel(), clon.ravel()
    area_acc = flow.accumulate(res["order"], rec, np.where(landf, A.ravel(), 0.0))
    rivers = [{"lat": round(float(LA[c]), 2), "lon": round(float(LO[c]), 2), "discharge_m3s": round(float(q[c])),
               "basin_km2": round(float(area_acc[c]))} for c in order_mouths[:25]]
    fid = np.load(args.canon_relief)["feature_id"].ravel()
    canon = {}
    for k, name in ((25, "Shrinking Sea"), (26, "Bleeding Sea"), (7, "The Neck"), (28, "Dry Deep")):
        m = landf & (fid == k)
        if m.any():
            canon[name] = {"lake_share": round(float(A.ravel()[m & np.isin(lake.ravel(), (1, 2))].sum() / A.ravel()[m].sum()), 3),
                           "marsh_share": round(float(A.ravel()[m & (lake.ravel() == 3)].sum() / A.ravel()[m].sum()), 3),
                           "salt_flat_share": round(float(A.ravel()[m & (lake.ravel() == 4)].sum() / A.ravel()[m].sum()), 3),
                           "max_discharge_m3s": round(float(q[m].max()))}
    law = landf & (LA > 70) & (LO > -35) & (LO < 25)
    canon["Land of Always Winter lakes"] = {"lake_km2": round(float(A.ravel()[law & np.isin(lake.ravel(), (1, 2))].sum()))}
    # the Rhoyne (E8): the largest river mouth within 6 degrees of Volantis
    vla, vlo = wc_places()["Volantis"]
    near = np.flatnonzero(mouth & (np.abs(LA - vla) < 4) & (np.abs(LO - vlo) < 6))
    if near.size:
        c = near[np.argmax(q[near])]
        canon["Rhoyne at Volantis"] = {"discharge_m3s": round(float(q[c])), "basin_km2": round(float(area_acc[c])),
                                       "lat": round(float(LA[c]), 2), "lon": round(float(LO[c]), 2)}
    # each canon river: discharge where it ends (largest discharge among its cells) and whether its cells drain
    # into each other (share of its cells whose receiver is a canon river cell or the sea)
    names = [r["name"] for r in json.loads((HERE.parent / "data" / "canon_rivers.json").read_text())["rivers"]]
    ridf = rid.ravel()
    canon_rivers = {}
    for k, nm in enumerate(names):
        cells = np.flatnonzero((ridf == k + 1) & landf)
        if cells.size == 0:
            continue
        r_ = rec[cells]
        follows = (r_ < 0) | (ridf[np.maximum(r_, 0)] > 0) | ~landf[np.maximum(r_, 0)]
        canon_rivers[nm] = {"max_discharge_m3s": round(float(q[cells].max())),
                            "basin_km2": round(float(area_acc[cells].max())),
                            "cells_draining_along_canon_rivers": round(float(follows.mean()), 3)}
    canon["rivers"] = canon_rivers
    canon["lakes"] = canon_lakes
    Al = A.ravel()
    stats = {
        "land_mean_runoff_mm": round(float((res["runoff_mm"] * Al)[landf].sum() / Al[landf].sum())),
        "runoff_to_sea_km3_yr": round(float(q[mouth].sum()) * flow.SECONDS_PER_YEAR / 1e9),
        "lake_area_km2": round(float(Al[lake.ravel() == 1].sum())),
        "terminal_lake_area_km2": round(float(Al[lake.ravel() == 2].sum())),
        "marsh_area_km2": round(float(Al[lake.ravel() == 3].sum())),
        "salt_flat_area_km2": round(float(Al[lake.ravel() == 4].sum())),
        "largest_lakes": sorted(res["lakes"], key=lambda x: -x["area_km2"])[:12],
        "largest_rivers_at_mouth": rivers,
        "canon": canon,
        "earth_reference": {"runoff_to_sea_km3_yr": 39000, "lakes_km2": 2.7e6},
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, H = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(H) + 0.5) * 180 / H
    PLO, PLA = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(PLA), np.radians(PLO), n)
    e = surface[ff, jj, ii]
    lnd = land[ff, jj, ii]
    shade = np.clip(np.gradient(e, axis=1) * 0.02 + 0.5, 0, 1)
    img = np.where(lnd[..., None], np.stack([170 + 50 * shade, 160 + 40 * shade, 130 + 30 * shade], -1),
                   np.array([40, 70, 120]))
    img = np.where((sheet[ff, jj, ii])[..., None], np.array([235, 240, 250]), img)
    qq = q.reshape(shape)[ff, jj, ii]
    lk = lake.reshape(shape)[ff, jj, ii]
    w = np.clip((np.log10(np.maximum(qq, 1)) - 1) / 3.5, 0, 1)
    riv = lnd & (qq > hp["channel_min_m3s"] * 3) & ~sheet[ff, jj, ii] & (lk == 0)
    img = np.where(riv[..., None], (1 - w[..., None]) * np.array([120, 170, 230]) + w[..., None] * np.array([10, 40, 170]), img)
    for code, col in ((3, (120, 160, 140)), (4, (235, 230, 215)), (1, (50, 100, 200)), (2, (90, 150, 200))):
        img = np.where((lk == code)[..., None], np.array(col), img)
    Image.fromarray(img.astype(np.uint8)).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


def carve_canon(surface, land, cr, hp, nb):
    """Carve the canon rivers so that they flow the way the map draws them. Each river, from its source, falls at
    least canon_drop_m per cell and lies at least canon_burn_m below the surface; it never goes below a floor of
    canon_drop_m x (cells to the sea along the river network), which keeps it above sea level and above the river
    it joins, so it cannot drain out through a low-lying source near the coast. A river that joins another ends in
    a cell of that river."""
    names = json.loads((HERE.parent / "data" / "canon_rivers.json").read_text())["rivers"]
    flat = surface.ravel().copy()
    base = surface.ravel() - hp["canon_burn_m"]
    landf = land.ravel()
    drop = hp["canon_drop_m"]

    def depth(r):
        d, j = 0, r.get("joins")
        while j:
            d += 1
            j = next((x.get("joins") for x in names if x["name"] == j), None)
        return d

    to_sea = np.full(flat.size, -1, np.int64)        # cells to the sea along the network
    for k in sorted(range(len(names)), key=lambda k: depth(names[k])):
        seq = list(cr[f"r{k}"])
        cut, tail = len(seq), 0
        for i, c in enumerate(seq):
            if not landf[c]:
                cut = i
                break
            if names[k].get("joins") and i > 0 and to_sea[c] >= 0:
                cut, tail = i + 1, int(to_sea[c])          # include the junction cell
                break
        body = seq[:cut]
        if not body:
            continue
        if not names[k].get("joins") and cut == len(seq):
            # the traced mouth may stop short of the model's coast: continue by the shortest way to the sea
            body = body + to_coast(body[-1], landf, nb, 40)
        m = len(body)
        floor = drop * (np.arange(m, 0, -1) - 1 + tail)
        h = np.inf
        for i, c in enumerate(body):
            h = min(h - drop, base[c])
            v = min(max(h, floor[i]), surface.ravel()[c])
            flat[c] = v if to_sea[c] < 0 else min(flat[c], v)
            if to_sea[c] < 0:
                to_sea[c] = m - 1 - i + tail
    return flat.reshape(surface.shape)


def to_coast(start, landf, nb, max_steps):
    """Land cells on the shortest 8-neighbour path from start to the nearest sea cell (excluding both ends)."""
    prev = {int(start): -1}
    frontier = [int(start)]
    for _ in range(max_steps):
        nxt = []
        for c in frontier:
            for q in nb[:, c]:
                q = int(q)
                if q in prev:
                    continue
                prev[q] = c
                if not landf[q]:
                    path = []
                    p = c
                    while p != int(start):
                        path.append(p)
                        p = prev[p]
                    return path[::-1]
                nxt.append(q)
        frontier = nxt
    return []


def wc_places():
    import climate as tlayer
    return tlayer.PLACES


if __name__ == "__main__":
    main()
