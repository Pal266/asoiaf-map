"""Layer "resources": mineral deposits and area resources (experiment E17).

The three-part method of USGS mineral resource assessments (Singer & Menzie 2010), run forward to make a world:
1. Permissive tracts. Each deposit type of the USGS descriptive models (Cox & Singer 1986) forms in a geological
   setting. Its tract is where the layers above show that setting: tectonic province, rock (lithology E14), age
   (geology E16), nearness to intrusions, faults or evaporites, and for laterites the climate. Placers lie on the
   alluvium of rivers within reach of lode-gold tracts.
2. Number of deposits. Each tract (permissive outcrops closer than 30 km joined, as tracts are delineated at
   regional map scale) gets the median number of deposits of Singer's (2008) general deposit
   density regression, log10(deposits per 100,000 km2) = 4.2474 - 0.5146 log10(area) - 0.2254 log10(median
   tonnage). The actual number is a Poisson draw around it.
3. Size and grade. Each deposit's tonnage is drawn from the type's lognormal grade-and-tonnage model, truncated at
   3 standard deviations (the models describe observed deposits, not an unbounded tail); grades are the model's
   medians.
Deposits are placed at random cells of their tract; vein types prefer cells within 20 km of a fault.

Area resources (no counts, the ground itself): coal measures, salt, bog iron, peat, clay, building stones, obsidian,
sulfur.

Input:
  --lithology, --geology, --temperature, --precipitation, --soils, --vegetation, --canon-relief (npz of those layers)
  world/resources/data/resources.json
Output (in --out):
  resources.npz [6, 1000, 1000]: deposit (int32, index into deposits.json, -1 none), permissive (uint32 bit per
                 deposit model, in the order of resources.json), area (uint16 bit per area resource)
  deposits.json  every deposit: type, cell, latitude and longitude, tonnage, grades, contained metal
  stats.json     counts and contained metal by type and region, canon checks
  preview.png
Deterministic.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
sys.path.insert(0, str(WORLD / "climate" / "scripts"))
sys.path.insert(0, str(WORLD / "hydrology" / "scripts"))
import cubesphere as cs     # noqa: E402
import world_climate as wc  # noqa: E402

REGIONS = {"Westeros south of the Wall": (28.0, 63.75, -30.0, 15.0), "Beyond the Wall": (63.75, 90.0, -40.0, 30.0),
           "Western Essos": (25.0, 55.0, 15.0, 50.0)}
IRON_ISLANDS = (45.0, 47.5, -16.0, -8.0)        # the islands of Ironman's Bay (lat, lat, lon, lon)


def dilate(mask, k, halo):
    F, J, I = halo
    m = mask.copy()
    for _ in range(int(k)):
        g = m[F, J, I]
        m = (g[:, :-2, 1:-1] | g[:, 2:, 1:-1] | g[:, 1:-1, :-2] | g[:, 1:-1, 2:] | g[:, 1:-1, 1:-1]
             | g[:, :-2, :-2] | g[:, :-2, 2:] | g[:, 2:, :-2] | g[:, 2:, 2:])
    return m


def components(mask, halo):
    """Connected components (8-neighbour, across face seams) of a cube mask: labels (-1 outside)."""
    F, J, I = halo
    n = mask.shape[1]
    idx = np.full(mask.shape, -1, np.int64)
    cells = np.flatnonzero(mask)
    idx.ravel()[cells] = np.arange(len(cells))
    g = idx[F, J, I]
    centre = g[:, 1:-1, 1:-1]
    rows, cols = [], []
    for a, b in ((0, 1), (1, 0), (1, 1), (1, -1)):
        nb = g[:, 1 + a:1 + a + n, 1 + b:1 + b + n]
        m = (centre >= 0) & (nb >= 0)
        rows.append(centre[m])
        cols.append(nb[m])
    r, c = np.concatenate(rows), np.concatenate(cols)
    k, lab = connected_components(coo_matrix((np.ones(len(r)), (r, c)), shape=(len(cells), len(cells))),
                                  directed=False)
    out = np.full(mask.shape, -1, np.int64)
    out.ravel()[cells] = lab
    return out, k


def main():
    ap = argparse.ArgumentParser()
    for a in ("--lithology", "--geology", "--temperature", "--precipitation", "--soils", "--vegetation",
              "--canon-relief", "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "resources.json").read_text())
    lcfg = json.loads((WORLD / "lithology" / "data" / "lithology.json").read_text())
    scfg = json.loads((WORLD / "soils" / "data" / "soils.json").read_text())
    vcfg = json.loads((WORLD / "vegetation" / "data" / "vegetation.json").read_text())
    n = cs.N
    halo = cs.halo_index(n, 1)
    L = np.load(args.lithology)
    G = np.load(args.geology)
    rock, deposit, prov = L["rock"], L["deposit"], L["province"]
    age = G["age_ma"].astype(np.float64)
    bed_age = age.copy()
    fault = (G["fault_type"] >= 3)                                   # transforms, thrusts, normal faults, sutures
    land = ~np.isin(rock, [lcfg["order"].index("wb"), lcfg["order"].index("ig")])
    T = np.load(args.temperature)["t_month_dC"].astype(np.float32) / 10
    P = np.load(args.precipitation)["p_month_mm"].astype(np.float32)
    mat, pyr = T.mean(0), P.sum(0)
    rsg = np.load(args.soils)["rsg"]
    biome = np.load(args.vegetation)["biome"]
    fid = np.load(args.canon_relief)["feature_id"]
    clat, clon = wc.cell_latlon_deg(n)
    A = np.broadcast_to(cs.cell_area_km2(n), land.shape)
    R = {c: i for i, c in enumerate(lcfg["order"])}

    def rocks(names):
        return np.isin(rock, [R[c] for c in names]) | np.isin(L["bedrock"], [R[c] for c in names])

    near_fault = {k: dilate(fault & land, k, halo) for k in (2, 3)}
    ftype = G["fault_type"]
    rng = np.random.default_rng(cfg["seed"])
    dens = cfg["density"]
    models = cfg["models"]
    permissive = np.zeros(land.shape, np.uint32)
    tracts = {}
    deposits = []
    dep_map = np.full(land.shape, -1, np.int32)
    order = list(models)
    for bit, key in enumerate(order):
        m = models[key]
        t = m["tract"]
        p = land.copy()
        if "province" in t:
            p &= np.isin(prov, t["province"])
        if "rock" in t:
            p &= rocks(t["rock"])
        if "deposit" in t:
            p &= np.isin(deposit, t["deposit"])
        if "age_max_ma" in t:
            p &= bed_age <= t["age_max_ma"]
        if "age_min_ma" in t:
            p &= bed_age >= t["age_min_ma"]
        if "near_rock" in t:
            p &= dilate(land & rocks(t["near_rock"]), t["near_cells"], halo)
        if "near_fault" in t:
            nf = t["near_fault"]
            p &= dilate(np.isin(ftype, nf["types"]) & land, nf["cells"], halo)
        if "near_province" in t:
            npv = t["near_province"]
            p &= dilate(np.isin(prov, npv["province"]) & land, npv["cells"], halo)
        if "climate" in t:
            p &= (mat >= t["climate"]["mat_min_c"]) & (pyr >= t["climate"]["p_min_mm"])
        if "near_types" in t:
            src = np.zeros(land.shape, bool)
            for k2 in t["near_types"]:
                src |= tracts.get(k2, np.zeros(land.shape, bool))
            p &= dilate(src, round(t["near_km"] / 10), halo)
        tracts[key] = p
        permissive |= (p.astype(np.uint32) << np.uint32(bit))
        # deposits, tract by tract: outcrops closer than tract_merge_km belong to one tract
        lab, k = components(dilate(p, round(cfg["tract_merge_km"] / 20), halo), halo)
        lab = np.where(p, lab, -1)
        if k == 0:
            continue
        flat = lab.ravel()
        cells = np.flatnonzero(flat >= 0)
        lab_c = flat[cells]
        area = np.bincount(lab_c, weights=A.ravel()[cells], minlength=k)
        logS = np.log10(m["tonnage_mt"])
        for tr in np.flatnonzero(area >= cfg["min_tract_km2"]):
            a = area[tr]
            n50 = a / 1e5 * 10 ** (dens["a"] + dens["b_area"] * np.log10(a) + dens["c_size"] * logS)
            nd = rng.poisson(n50)
            if nd == 0:
                continue
            tc = cells[lab_c == tr]
            w = np.ones(len(tc))
            if m["place"] == "faults":
                w += 4 * near_fault[2].ravel()[tc]
            nd = min(nd, len(tc))
            chosen = rng.choice(tc, size=nd, replace=False, p=w / w.sum())
            for c in chosen:
                f, j, i = np.unravel_index(c, land.shape)
                z = np.clip(rng.standard_normal(), -cfg["tonnage_z_max"], cfg["tonnage_z_max"])
                ton = 10 ** (logS + m["sd_log10"] * z)
                metal = {}
                for el, g in m["commodities"].items():
                    unit = m["units"][el]
                    metal[el] = ton * 1e6 * (g / 100 if unit == "%" else g / 1e6)     # tonnes of metal
                if dep_map[f, j, i] < 0:
                    dep_map[f, j, i] = len(deposits)
                deposits.append({"type": key, "cell": [int(f), int(j), int(i)],
                                 "lat": round(float(clat[f, j, i]), 3), "lon": round(float(clon[f, j, i]), 3),
                                 "tonnage_mt": round(float(ton), 4),
                                 "grades": m["commodities"], "units": m["units"],
                                 "contained_t": {el: round(float(v), 2) for el, v in metal.items()}})

    # ---- area resources ----
    ar = cfg["area_resources"]
    G2 = {g: i for i, g in enumerate(scfg["groups"])}
    B = {b: i for i, b in enumerate(vcfg["biomes"])}
    area_bits = np.zeros(land.shape, np.uint16)
    area_names = [k for k in ar if not k.startswith("_")]
    for bit, key in enumerate(area_names):
        r = ar[key]
        mm = np.zeros(land.shape, bool)
        if "rock" in r:
            mm |= rocks(r["rock"]) & land
        if "deposit" in r:
            mm |= np.isin(deposit, r["deposit"])
        if "rock_or" in r:
            mm |= rocks(r["rock_or"]) & land
        if "soils" in r:
            mm |= np.isin(rsg, [G2[s] for s in r["soils"]])
        if "rock" in r and "deposit" in r:
            mm = (rocks(r["rock"]) & land) | np.isin(deposit, r["deposit"])
        if "age_windows_ma" in r:
            mm &= np.any([(bed_age >= lo) & (bed_age <= hi) for lo, hi in r["age_windows_ma"]], axis=0)
        if "age_max_ma" in r:
            mm &= bed_age <= r["age_max_ma"]
        if "province" in r:
            mm &= np.isin(prov, r["province"])
        if "near_rock" in r:
            mm &= dilate(land & rocks(r["near_rock"]), r["near_cells"], halo)
        if "vegetation_not" in r:
            mm &= ~np.isin(biome, [B[b] for b in r["vegetation_not"]])
        area_bits |= (mm.astype(np.uint16) << np.uint16(bit))
    np.savez_compressed(out / "resources.npz", deposit=dep_map, permissive=permissive, area=area_bits)
    (out / "deposits.json").write_text(json.dumps({"models": order, "area_resources": area_names,
                                                   "deposits": deposits}, ensure_ascii=False))

    # ---- statistics ----
    def summary(sel):
        s = {}
        for d in sel:
            e = s.setdefault(d["type"], {"count": 0, "tonnage_mt": 0.0, "contained_t": {}})
            e["count"] += 1
            e["tonnage_mt"] += d["tonnage_mt"]
            for el, v in d["contained_t"].items():
                e["contained_t"][el] = e["contained_t"].get(el, 0) + v
        for e in s.values():
            e["tonnage_mt"] = round(e["tonnage_mt"], 2)
            e["contained_t"] = {el: round(v) for el, v in e["contained_t"].items()}
        return s

    def metal_totals(sel):
        tot = {}
        for d in sel:
            for el, v in d["contained_t"].items():
                tot[el] = tot.get(el, 0) + v
        return {el: round(v) for el, v in sorted(tot.items())}

    def inside(d, box):
        a0, a1, o0, o1 = box
        return a0 <= d["lat"] < a1 and o0 <= d["lon"] < o1

    regions = {name: {"metals_t": metal_totals([d for d in deposits if inside(d, box)]),
                      "counts": {k: v["count"] for k, v in summary([d for d in deposits if inside(d, box)]).items()}}
               for name, box in REGIONS.items()}
    lnd = land
    permissive_km2 = {k: round(float(A[tracts[k]].sum())) for k in order}
    area_km2 = {k: round(float(A[(area_bits >> np.uint16(b)) & 1 == 1].sum())) for b, k in enumerate(area_names)}

    def canon_box(mask):
        sel = [d for d in deposits if mask[tuple(d["cell"])]]
        return {"deposits": {k: v["count"] for k, v in summary(sel).items()}, "metals_t": metal_totals(sel),
                "area_resources_share": {k: round(float(A[mask & ((area_bits >> np.uint16(b)) & 1 == 1)].sum()
                                                        / A[mask].sum()), 3)
                                         for b, k in enumerate(area_names)}}

    def point_mask(f, radius_km):
        """Land within radius_km of a canon feature drawn as a point."""
        sys.path.insert(0, str(WORLD / "canon_outline" / "scripts"))
        import canon_outline as co
        la, lo = co.globe(*f["point_px"])
        dlat, dlon = np.radians(clat) - la, np.radians(clon) - lo
        hav = np.sin(dlat / 2) ** 2 + np.cos(la) * np.cos(np.radians(clat)) * np.sin(dlon / 2) ** 2
        return 2 * cs.R_KM * np.arcsin(np.sqrt(np.minimum(hav, 1))) < radius_km

    a0, a1, o0, o1 = IRON_ISLANDS
    ii = lnd & (clat >= a0) & (clat < a1) & (clon >= o0) & (clon < o1)
    feats = json.loads((WORLD / "canon_relief" / "data" / "canon_features.json").read_text())["features"]
    fnum = {f["name"]: k + 1 for k, f in enumerate(feats)}
    west = lnd & (fid == fnum["Westerlands hills (Golden Tooth)"])
    west_near = dilate(west, 10, halo) & lnd
    canon = {"Westerlands hills (and 100 km around)": canon_box(west_near),
             "Iron Islands": canon_box(ii),
             "Fourteen Flames and the Valyrian peninsula": canon_box(lnd & dilate(
                 fid == fnum["Fourteen Flames and the Valyrian peninsula"], 5, halo)),
             "Dragonstone (30 km around Dragonmont)": canon_box(lnd & point_mask(feats[fnum["Dragonmont (Dragonstone)"] - 1],
                                                                                  30.0))}
    stats = {"deposits_total": len(deposits), "by_type": summary(deposits), "permissive_km2": permissive_km2,
             "area_resources_km2": area_km2, "regions": regions, "canon": canon,
             "earth_known_counts": cfg["earth_known_counts"]}
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, H = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(H) + 0.5) * 180 / H
    PLO, PLA = np.meshgrid(plon, plat)
    ff, jj, i2 = cs.latlon_to_cell(np.radians(PLA), np.radians(PLO), n)
    img = np.where(land[ff, jj, i2][..., None], np.array([225, 220, 200], np.uint8), np.array([70, 120, 195], np.uint8))
    img = img.astype(np.uint8)
    img[((area_bits >> np.uint16(area_names.index("coal"))) & 1 == 1)[ff, jj, i2]] = (90, 90, 90)
    col = {"Au": (230, 180, 0), "Ag": (180, 180, 220), "Cu": (200, 100, 40), "Fe": (150, 40, 40), "Sn": (60, 160, 160),
           "Zn": (120, 60, 160), "Pb": (60, 60, 120), "Ni": (0, 150, 0), "Hg": (220, 0, 120)}
    for d in deposits:
        el = next(iter(d["grades"]))
        x = int((d["lon"] + 180) / 360 * W) % W
        y = min(max(int((90 - d["lat"]) / 180 * H), 0), H - 1)
        img[max(y - 1, 0):y + 2, max(x - 1, 0):x + 2] = col.get(el, (0, 0, 0))
    Image.fromarray(img).save(out / "preview.png")
    print(json.dumps({k: stats[k] for k in ("deposits_total", "permissive_km2", "area_resources_km2", "canon")},
                     indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
