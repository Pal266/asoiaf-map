"""Layer "ice": ice sheets and glaciers of the world in equilibrium with its climate (experiment E11).

The ice and the climate are found together, by iteration from bare ground:
1. climate of the current surface: the temperature model (E09) with the land height and the ice sheets' albedo,
   and the moisture model (E10) - world/climate/scripts/world_climate.py;
2. surface mass balance of every 10 km cell (positive degree days, smb.py) at the ice surface;
3. ice sheets: connected areas of positive balance, at least sheet_min_km2, grown into the neighbouring ablation
   zone, least negative balance first, for as long as the sheet's total balance stays positive (an equilibrium
   sheet's accumulation equals its ablation; no calving: the sheets stop at the coast);
4. ice thickness: a perfectly plastic ice sheet (Nye 1951; the basal shear stress tau0 everywhere at its yield
   value): H |grad s| = tau0 / (rho g), integrated inwards from the margin over the real bed (Dijkstra order) with
   local isostatic depression of the bed by rho_ice / rho_mantle of the thickness;
5. mountain glaciers outside the sheets: the share of each cell above the equilibrium line, from its sub-grid
   relief (smb.glacier_fraction).
Repeat with the new surface until the sheets stop changing.

Input:
  --relief   relief.npz (elevation_m, subgrid_relief_m) from world/relief/scripts/relief.py
  world/ice/data/ice.json, world/climate/data/climate.json, world/precipitation/data/precipitation.json
Output (in --out):
  ice.npz     [6, 1000, 1000]: ice_thickness_m float32, surface_m float32 (ice surface on land, bedrock where no
              ice, after isostasy; sea floor at sea), sheet_frac uint8 (255 = ice sheet; else glacier share x 255)
  stats.json  areas, volumes, iterations, canon checks
  preview.png ice thickness over shaded relief
Deterministic.
"""

import argparse
import heapq
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
import smb as smbm          # noqa: E402

G = 9.81
CANON = {"Land of Always Winter": (70.0, 90.0, -35.0, 25.0)}   # rough box north of the Frostfangs (lat, lon)


def neighbours(n=cs.N):
    """Flat indices of the 4 edge neighbours of every cell (4, 6 n n), across face seams."""
    F, J, I = cs.halo_index(n, 1)
    flat = (F * n * n + J * n + I).astype(np.int32)
    return np.stack([flat[:, 1:-1, 2:], flat[:, 1:-1, :-2], flat[:, 2:, 1:-1], flat[:, :-2, 1:-1]]).reshape(4, -1)


def components(mask, nb):
    """Connected components of a flat bool mask (4-neighbour). Returns a label array (-1 outside)."""
    lab = np.full(mask.size, -1, np.int32)
    k = 0
    for s in np.flatnonzero(mask):
        if lab[s] >= 0:
            continue
        stack = [s]
        lab[s] = k
        while stack:
            c = stack.pop()
            for q in nb[:, c]:
                if mask[q] and lab[q] < 0:
                    lab[q] = k
                    stack.append(q)
        k += 1
    return lab, k


def calving_loss(land, sea_depth, step, nb, cfg):
    """Mass a land cell would lose each year by calving if an ice sheet reached the sea there (m w.e. x km2):
    calving speed proportional to water depth, u = c D (Brown, Meier & Post 1982), through the plastic margin
    thickness, across every edge the cell shares with the sea."""
    sh = cfg["sheet"]
    r = sh["rho_ice"] / sh["rho_mantle"]
    c = (1 - r) * sh["tau0_kPa"] * 1e3 / (sh["rho_ice"] * G)
    land = land.ravel()
    loss = np.zeros(land.size)
    for k in range(nb.shape[0]):
        q = nb[k]
        edge = land & ~land[q]
        D = np.minimum(sea_depth.ravel()[q], sh["calving_depth_max_m"])
        Hm = np.sqrt(c * step[k]) / (1 - r)
        width = step[2 if k < 2 else 0]          # edge length: the step in the other direction
        loss += np.where(edge, sh["rho_ice"] / 1000 * sh["calving_c_per_yr"] * D * Hm * width, 0) / 1e6
    return loss


def sheets(b, area, land, nb, min_km2, calv):
    """Ice-sheet mask: components of b > 0 of at least min_km2 whose balance minus calving is positive, grown while
    their total balance (surface balance minus calving where they reach the sea) stays >= 0."""
    lab, k = components((b > 0) & land, nb)
    sheet = np.zeros(b.size, bool)
    for c in range(k):
        cells = np.flatnonzero(lab == c)
        total = float((b[cells] * area[cells] - calv[cells]).sum())
        if area[cells].sum() < min_km2 or total < 0:
            continue
        sheet[cells] = True
        heap = []
        seen = set(cells.tolist())
        for cc in cells:
            for q in nb[:, cc]:
                if land[q] and not sheet[q] and q not in seen:
                    seen.add(q)
                    heapq.heappush(heap, (-(b[q] * area[q] - calv[q]) / area[q], q))
        while heap:
            nbq, q = heapq.heappop(heap)
            gain = b[q] * area[q] - calv[q]
            if total + gain < 0:
                break
            total += gain
            sheet[q] = True
            for r in nb[:, q]:
                if land[r] and not sheet[r] and r not in seen:
                    seen.add(r)
                    heapq.heappush(heap, (-(b[r] * area[r] - calv[r]) / area[r], r))
    return sheet


def plastic_thickness(sheet, bed, step, nb, tau0, rho_i, rho_m):
    """Perfectly plastic ice over the bed, with local isostasy. step (4, cells): distance (m) from each cell to
    each of its neighbours. Returns thickness H (m) on the flat grid."""
    r = rho_i / rho_m
    c = (1 - r) * tau0 / (rho_i * G)          # h |grad s| = c, h = (1 - r) H the height above the unloaded bed
    s = np.full(bed.size, np.inf)
    heap = []
    for q in np.flatnonzero(sheet):
        out = ~sheet[nb[:, q]]
        if out.any():                         # margin cell: half a step from the edge
            s[q] = bed[q] + np.sqrt(c * step[out, q].min())
            heapq.heappush(heap, (s[q], q))
    done = np.zeros(bed.size, bool)
    while heap:
        sq, q = heapq.heappop(heap)
        if done[q]:
            continue
        done[q] = True
        for k in range(nb.shape[0]):
            p = nb[k, q]
            if not sheet[p] or done[p]:
                continue
            cand = bed[p] + np.sqrt(max(sq - bed[p], 0.0) ** 2 + 2 * c * step[k, q])
            if cand < s[p]:
                s[p] = cand
                heapq.heappush(heap, (cand, p))
    h = np.where(sheet, np.maximum(s - bed, 0), 0)
    return h / (1 - r)


def equilibrate(geom, cfg, start="bare", verbose=True):
    """Find ice and climate together on any grid. geom: dict with
    land, bed, relief, A (km2) - arrays of the grid's shape; nb (4, cells) flat neighbour indices; step (4, cells)
    distance (m) to each neighbour; sea_depth (m, 0 on land);
    climate(surface, albedo_share or None) -> (T (12, ...) C at the surface, P (12, ...) mm)."""
    land, bed, relief, A, step, nb = (geom[k] for k in ("land", "bed", "relief", "A", "step", "nb"))
    p, g, sh = cfg["smb"], cfg["glacier"], cfg["sheet"]
    r = sh["rho_ice"] / sh["rho_mantle"]
    tau = sh["tau0_kPa"] * 1e3

    def thickness(mask):
        return plastic_thickness(mask.ravel(), bed.ravel(), step, nb, tau,
                                 sh["rho_ice"], sh["rho_mantle"]).reshape(land.shape)

    calv = calving_loss(land, geom["sea_depth"], step, nb, cfg)
    H = np.zeros(land.shape)
    sheet = np.zeros(land.shape, bool)
    frac = np.zeros(land.shape)
    if start in ("glacial", "history"):
        # start from a glacial state: ice sheets wherever the bare-ground annual mean temperature is below
        # glacial_start_T_C (a cold-period ice cover), then let the present climate keep or remove them.
        # "history": only inside the regions the geological history has kept polar (geom["glacial_regions"])
        T0, _ = geom["climate"](np.where(land, bed, 0.0), None)
        sheet = land & (T0.mean(0) < sh["glacial_start_T_C"])
        if start == "history":
            sheet &= geom["glacial_regions"]
        H = thickness(sheet)
    history = []
    for it in range(sh["iterations"]):
        surface = np.where(land, bed + (1 - r) * H, 0.0)
        albedo_share = np.maximum(sheet.astype(float), frac)
        T, P = geom["climate"](surface, albedo_share if albedo_share.any() else None)
        # inner loop: with this climate, the ice and its own height (lapse rate, elevation desertification) are
        # brought to equilibrium; the climate's large-scale response (albedo, circulation) waits for the next pass
        Hn, new_sheet = H, sheet
        for _ in range(sh["inner_iterations"]):
            dz = np.where(land, bed + (1 - r) * Hn - surface, 0.0)
            b0, parts = smbm.smb(T, P, p, dz=dz)
            bu, _ = smbm.smb(T, P, p, dz=dz + g["dz_up_m"])
            ns = sheets(b0.ravel(), A.ravel(), land.ravel(), nb, sh["sheet_min_km2"], calv).reshape(land.shape)
            Hn = thickness(ns)
            inner_change = float(A[ns ^ new_sheet].sum())
            new_sheet = ns
            if inner_change < sh["converged_km2"]:
                break
        frac = np.where(land & ~new_sheet, smbm.glacier_fraction(b0, bu, g["dz_up_m"], relief, g["aar"]), 0.0)
        changed = float(A[new_sheet ^ sheet].sum())
        sheet, H = new_sheet, Hn
        rec = {"iteration": it, "sheet_Mkm2": round(float(A[sheet].sum()) / 1e6, 3),
               "glacier_Mkm2": round(float((A * frac).sum()) / 1e6, 3), "changed_km2": round(changed),
               "volume_Mkm3": round(float((A * H).sum()) / 1e9, 3)}
        history.append(rec)
        if verbose:
            print(rec, flush=True)
        if it > 0 and changed < sh["outer_converged_km2"]:
            break
    return {"H": H, "surface": np.where(land, bed + (1 - r) * H, 0.0), "sheet": sheet, "frac": frac, "b0": b0,
            "T": T, "P": P, "history": history}


def run(elev, relief, cfg, start="bare", verbose=True):
    n = cs.N
    land = elev > 0
    clat, clon = wc.cell_latlon_deg(n)
    A = np.broadcast_to(cs.cell_area_km2(n), land.shape)

    def climate(surface, share):
        temp = wc.temperature(surface, land, share, clat, clon)
        prec = wc.precipitation(temp, surface, land, clat, clon)
        return temp["t_cube"], prec["p_cube"]

    nb = neighbours(n)
    dx = (np.sqrt(A) * 1000).ravel()
    geom = {"land": land, "bed": np.where(land, elev, 0.0), "relief": relief, "A": A,
            "sea_depth": np.where(land, 0.0, -elev),
            "step": 0.5 * (dx[None, :] + dx[nb]), "nb": nb, "climate": climate}
    reg = np.zeros(land.shape, bool)
    for r in cfg["sheet"]["history_glacial_regions"]:
        reg |= (clat >= r["lat_min"]) & (clat <= r["lat_max"])
    geom["glacial_regions"] = reg
    res = equilibrate(geom, cfg, start, verbose)
    res["surface"] = np.where(land, res["surface"], elev)
    res.update(clat=clat, clon=clon, A=A, land=land)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--relief", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", choices=["bare", "glacial", "history"], default=None,
                    help="initial state (default: ice.json sheet.start)")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "ice.json").read_text())
    rel = np.load(args.relief)
    elev = rel["elevation_m"].astype(np.float64)
    start = args.start or cfg["sheet"]["start"]
    res = run(elev, rel["subgrid_relief_m"].astype(np.float64), cfg, start)
    H, sheet, frac, land, A = res["H"], res["sheet"], res["frac"], res["land"], res["A"]
    sf = np.where(sheet, 255, np.round(frac * 254)).astype(np.uint8)
    np.savez_compressed(out / "ice.npz", ice_thickness_m=H.astype(np.float32),
                        surface_m=res["surface"].astype(np.float32), sheet_frac=sf)

    clat, clon = res["clat"], res["clon"]
    canon = {}
    for name, (a0, a1, o0, o1) in CANON.items():
        m = land & (clat >= a0) & (clat <= a1) & (clon >= o0) & (clon <= o1)
        canon[name] = {"land_km2": round(float(A[m].sum())),
                       "ice_sheet_share": round(float(A[m & sheet].sum() / max(A[m].sum(), 1)), 3),
                       "glacier_share_outside_sheets": round(float((A * frac)[m].sum() / max(A[m].sum(), 1)), 3),
                       "max_thickness_m": int(H[m].max()) if m.any() else 0}
    stats = {
        "start": start,
        "iterations": res["history"],
        "ice_sheet_area_Mkm2": round(float(A[sheet].sum()) / 1e6, 3),
        "glacier_area_Mkm2": round(float((A * frac).sum()) / 1e6, 3),
        "ice_volume_Mkm3": round(float((A * H).sum()) / 1e9, 3),
        "max_thickness_m": int(H.max()),
        "sea_level_equivalent_m": round(float((A * H).sum()) * 1e6 * 0.91 / 3.61e14, 1),
        "canon": canon,
        "earth_reference": {"greenland_Mkm2": 1.71, "greenland_max_thickness_m": 3200,
                            "antarctica_Mkm2": 14.0, "antarctica_max_thickness_m": 4800},
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, Hh = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(Hh) + 0.5) * 180 / Hh
    LON, LAT = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(LAT), np.radians(LON), cs.N)
    e = res["surface"][ff, jj, ii]
    lnd = land[ff, jj, ii]
    shade = np.clip(np.gradient(e, axis=1) * 0.02 + 0.5, 0, 1)
    img = np.where(lnd[..., None], (np.stack([150 + 60 * shade, 130 + 50 * shade, 100 + 40 * shade], -1)),
                   np.array([40, 70, 120]))
    h = H[ff, jj, ii]
    gl = frac[ff, jj, ii]
    ice_c = np.stack([np.interp(h, [0, 1000, 3000], c) for c in ([225, 170, 120], [240, 205, 170], [255, 245, 235])], -1)
    img = np.where((h > 0)[..., None], ice_c, img)
    img = np.where(((h <= 0) & (gl > 0.2))[..., None], np.array([210, 240, 255]), img)
    Image.fromarray(img.astype(np.uint8)).save(out / "preview.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
