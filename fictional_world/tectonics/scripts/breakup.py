"""E02: the plates at the breakup, 100 Ma (fictional_world/tectonics/METHOD.md, step 2).

  * Continental plates: the supercontinent splits along assembly sutures between blocks. Shorter sutures reactivate
    first (Erdős et al. 2025); the longest contacts stay welded until n_continental_plates remain.
  * Oceanic plates: the exterior ocean is divided among n_oceanic_plates (warped Voronoi cells of spaced seeds).
  * Euler poles: each continental plate moves away from the supercontinent's centre (dispersal), faster further out
    and slower with a larger cratonic share (Zahirovic et al. 2015); each oceanic plate moves towards the
    supercontinent (the exterior ocean is consumed at its margins). Every plate also gets a small random spin.
  * Floor of the exterior ocean: age from the distance to the nearest divergent boundary between oceanic plates and
    that boundary's half spreading rate, capped at max_ocean_age_myr.

Usage: python3 breakup.py WORK_DIR   (reads supercontinent.npz written by supercontinent.py)
"""

import os
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import supercontinent as sc  # noqa: E402
from supercontinent import CFG, Grid, cs, noise  # noqa: E402

CM_PER_YR = 0.1   # 1 km/Myr = 0.1 cm/yr


def tangent(c, target):
    """Unit tangent at c pointing along the great circle towards target."""
    d = target - (target @ c) * c
    n = np.linalg.norm(d)
    return d / n if n > 1e-9 else None


def pole_for(c, d, speed_cm, spin_deg, rng):
    """Rotation vector (rad/Myr) giving velocity speed_cm along d at c, plus a spin about c."""
    w = np.cross(c, d) * (speed_cm / CM_PER_YR / cs.R_KM)
    return w + c * np.radians(rng.normal(0, spin_deg))


def velocity(w, xyz):
    """cm/yr, shape (..., 3)."""
    return np.cross(w, xyz) * cs.R_KM * CM_PER_YR


def merge_blocks(grid, s, bk):
    """Continental plates as unions of blocks: weld the longest block contacts until the target number remains."""
    ter = s["terrane"]
    blk = np.where(ter >= 0, s["terrane_block"][np.maximum(ter, 0)], -1)
    n_b = blk.max() + 1
    cont = sc.contacts(grid, blk)
    parent = list(range(n_b))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    groups = n_b
    for (p, q), _ in sorted(cont.items(), key=lambda kv: -kv[1]):
        if groups <= bk["n_continental_plates"]:
            break
        rp, rq = find(p), find(q)
        if rp != rq:
            parent[rp] = rq
            groups -= 1
    roots = sorted({find(b) for b in range(n_b)})
    block_plate = np.array([roots.index(find(b)) for b in range(n_b)])
    rifted = [(p, q, c) for (p, q), c in cont.items() if block_plate[p] != block_plate[q]]
    welded = [(p, q, c) for (p, q), c in cont.items() if block_plate[p] == block_plate[q]]
    return blk, block_plate, rifted, welded


def microcontinents(grid, s, plate, n_c, bk, rng):
    """Small accreted terranes on a rifted edge break away as their own plates, as Madagascar, the Seychelles or the
    Lord Howe Rise left Pangaea's fragments. Candidates: non-cratonic terranes smaller than micro_max_area_km2 that
    touch a terrane of another continental plate, at most one per plate, never the last terrane of a plate."""
    ter = s["terrane"]
    area = grid.area
    n_t = ter.max() + 1
    tarea = np.bincount(ter[ter >= 0], weights=area[ter >= 0], minlength=n_t)
    tplate = np.array([np.bincount(plate[ter == t]).argmax() for t in range(n_t)])
    tc = sc.contacts(grid, ter)
    rift_edge = np.zeros(n_t, bool)
    for (a, b) in tc:
        if tplate[a] != tplate[b]:
            rift_edge[a] = rift_edge[b] = True
    cand = [t for t in range(n_t) if rift_edge[t] and not s["terrane_craton"][t]
            and tarea[t] < bk["micro_max_area_km2"] and (tplate == tplate[t]).sum() > 1]
    rng.shuffle(cand)
    used, micro = set(), []
    for t in cand:
        if len(micro) >= bk["n_microcontinents"] * bk["micro_candidates_factor"]:
            break
        if tplate[t] in used:
            continue
        used.add(tplate[t])
        plate = np.where(ter == t, n_c, plate)
        micro.append((int(t), int(n_c), int(tplate[t]), round(float(tarea[t]) / 1e6, 2)))
        n_c += 1
    return plate, n_c, micro


def oceanic_plates(grid, rng, ocean, bk, first_id):
    w = (grid.area * ocean).reshape(-1)
    pts = grid.xyz.reshape(-1, 3)
    seeds = []
    mind = np.cos(np.radians(bk["oceanic_seed_min_distance_deg"]))
    while len(seeds) < bk["n_oceanic_plates"]:
        p = pts[rng.choice(w.size, p=w / w.sum())]
        if all(p @ q < mind for q in seeds):
            seeds.append(p)
    seeds = np.array(seeds)
    warp = np.stack([noise(grid, rng, 25, 2, 2) for _ in range(3)], -1) * np.radians(6.0)
    p = grid.xyz + warp
    p /= np.linalg.norm(p, axis=-1, keepdims=True)
    lab = (p.reshape(-1, 3) @ seeds.T).argmax(1).reshape(ocean.shape) + first_id
    # keep only the largest piece of each plate; strays join a neighbour
    lab = np.where(ocean, lab, -1)
    for t in range(first_id, first_id + len(seeds)):
        comp = grid.components(lab == t)
        ids, cnt = np.unique(comp[comp >= 0], return_counts=True)
        for q in ids[cnt < cnt.max()] if len(ids) > 1 else []:
            cells = np.flatnonzero(comp == q)
            nb = lab.reshape(-1)[grid.nb4[cells].reshape(-1)]
            nb = nb[(nb >= first_id) & (nb != t)]
            if nb.size:
                lab.reshape(-1)[cells] = np.bincount(nb).argmax()
    return lab


def euler_poles(grid, plate, cont, crat, bk, rng):
    """Rotation vectors (rad/Myr) of all plates.

    Continental crust disperses: a plate with continental crust moves away from the continental crust of all other
    plates (inverse-square repulsion summed over their cells), which after a breakup points away from the
    supercontinent and later spreads the continents over the sphere. Oceanic plates move towards the continental
    crust of other plates, where they are subducted. Speed: oceanic plates oceanic_speed_cm_yr; plates with
    continental crust continental_speed_cm_yr, scaled by the strength of the repulsion and slowed by their cratonic
    share (Zahirovic et al. 2015), blended towards the oceanic speed by their oceanic share. Plus a random spin."""
    xyz = grid.xyz.reshape(-1, 3)
    area = grid.area.reshape(-1)
    pl = plate.reshape(-1)
    co = cont.reshape(-1)
    cr = crat.reshape(-1)
    n_p = pl.max() + 1
    cells = np.flatnonzero(co)
    W = np.zeros((n_p, 3))
    rep = {}
    cents = {}
    for p in np.unique(pl):
        m = pl == p
        c = (xyz[m] * area[m][:, None]).sum(0)
        c /= np.linalg.norm(c)
        cents[p] = c
        o = cells[pl[cells] != p]
        x = xyz[o]
        ang = np.arccos(np.clip(x @ c, -1, 1))
        t = x - (x @ c)[:, None] * c
        t /= np.maximum(np.linalg.norm(t, axis=1, keepdims=True), 1e-12)
        rep[p] = -(t * (area[o] / np.maximum(ang, np.radians(5)) ** 2)[:, None]).sum(0)
    cont_plates = [p for p in rep if area[(pl == p) & co].sum() > 0.25 * area[pl == p].sum()]
    ref = np.median([np.linalg.norm(rep[p]) for p in cont_plates]) if cont_plates else 1.0
    info = []
    for p, c in cents.items():
        m = pl == p
        fc = area[m & co].sum() / area[m].sum()
        frac = area[m & cr].sum() / max(area[m & co].sum(), 1e-9)
        if p in cont_plates:
            r = rep[p]
            lo, hi = bk["continental_speed_cm_yr"]
            spd = (lo + (hi - lo) * min(1, np.linalg.norm(r) / ref)) * (1 - bk["craton_slowdown"] * frac)
            spd = fc * spd + (1 - fc) * np.mean(bk["oceanic_speed_cm_yr"])
            kind = "continental"
        else:
            r = -rep[p]
            spd = rng.uniform(*bk["oceanic_speed_cm_yr"])
            kind = "oceanic"
            frac = 0.0
        d = tangent(c, c + r / max(np.linalg.norm(r), 1e-12) * 1e-3)
        if d is None:
            d = tangent(c, rng.standard_normal(3))
        W[p] = pole_for(c, d, spd, bk["spin_sd_deg_per_myr"], rng)
        info.append((p, kind, c, fc, frac))
    return W, info


def boundaries(grid, plate, W, thr):
    """Boundary edges between plates: cell, neighbour, normal and tangential relative speed (cm/yr) of the
    neighbour's plate with respect to the cell's plate; positive normal = diverging."""
    a = plate.reshape(-1)
    xyz = grid.xyz.reshape(-1, 3)
    nb = a[grid.nb4]
    cell, k = np.nonzero(nb != a[:, None])
    other = grid.nb4[cell, k]
    pa, pb = a[cell], a[other]
    r = xyz[cell]
    n = xyz[other] - xyz[cell]
    n -= (n * r).sum(1, keepdims=True) * r
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    rel = velocity(W[pb], r) - velocity(W[pa], r)
    vn = (rel * n).sum(1)
    vt = np.linalg.norm(rel - vn[:, None] * n, axis=1)
    kind = np.where(vn > thr, 1, np.where(vn < -thr, -1, 0))   # 1 divergent, -1 convergent, 0 transform
    kind = np.where((kind != 0) & (np.abs(vn) < 0.5 * vt), 0, kind)
    return dict(cell=cell, other=other, pa=pa, pb=pb, vn=vn, vt=vt, kind=kind)


def main(work):
    s = dict(np.load(os.path.join(work, "supercontinent.npz")))
    bk = CFG["breakup"]
    rng = np.random.default_rng(CFG["seed"] + 2)
    grid = Grid(CFG["sim_n"])
    cont = s["continent"]
    xyz = grid.xyz

    centre = (xyz[cont] * grid.area[cont][:, None]).sum(0)
    centre /= np.linalg.norm(centre)

    blk, block_plate, rifted, welded = merge_blocks(grid, s, bk)
    n_c = block_plate.max() + 1
    plate = np.where(cont, block_plate[np.maximum(blk, 0)], -1)
    plate, n_c, micro = microcontinents(grid, s, plate, n_c, bk, rng)
    plate = np.where(cont, plate, oceanic_plates(grid, rng, ~cont, bk, n_c))

    ter = s["terrane"]
    crat = np.where(ter >= 0, s["terrane_craton"][np.maximum(ter, 0)], False)
    W, info = euler_poles(grid, plate, cont, crat, bk, rng)
    # keep the microcontinents that move away from the plate they left (a sliver pushed back would weld again at
    # once), the fastest-opening first; the others stay with their plate, and the poles are recomputed
    if micro:
        b0 = boundaries(grid, plate, W, bk["boundary_threshold_cm_yr"])
        opening = []
        for t, mp, pp, a in micro:
            m = ((b0["pa"] == mp) & (b0["pb"] == pp)) | ((b0["pa"] == pp) & (b0["pb"] == mp))
            opening.append(b0["vn"][m].mean() if m.any() else -np.inf)
        order = np.argsort(opening)[::-1]
        keep = [micro[i] for i in order[:bk["n_microcontinents"]] if opening[i] > 0]
        for t, mp, pp, a in micro:
            if (t, mp, pp, a) not in keep:
                plate[plate == mp] = pp
        print("microcontinents (terrane, plate, parent plate, M km2):", keep,
              "; opening cm/yr:", [round(float(opening[i]), 1) for i in order])
        W, info = euler_poles(grid, plate, cont, crat, bk, rng)

    # RMS speeds
    total = grid.area.sum()
    print(f"plates: {sum(x[1] == 'continental' for x in info)} continental, {sum(x[1] == 'oceanic' for x in info)} oceanic")
    rms_all = []
    for p, kind, c, fc, frac in info:
        m = plate == p
        v = np.linalg.norm(velocity(W[p], xyz[m]), axis=1)
        rms = np.sqrt((grid.area[m] * v ** 2).sum() / grid.area[m].sum())
        rms_all.append((rms, grid.area[m].sum()))
        pole = W[p] / np.linalg.norm(W[p])
        plat, plon = np.degrees(cs.xyz_to_latlon(*pole))
        print(f"  plate {p:2d} {kind:11s} area {grid.area[m].sum() / total:5.3f}, craton {frac:4.2f}, "
              f"continental {fc:4.2f}, pole ({plat:6.1f}, {plon:7.1f}), "
              f"{np.degrees(np.linalg.norm(W[p])):5.3f} deg/Myr, RMS {rms:4.1f} cm/yr")
    r, a = np.array(rms_all).T
    print(f"area-weighted RMS speed {np.sqrt((a * r ** 2).sum() / a.sum()):.1f} cm/yr; continental plates median "
          f"{np.median(r[[i for i, x in enumerate(info) if x[1] == 'continental']]):.1f}, oceanic median {np.median(r[[i for i, x in enumerate(info) if x[1] == 'oceanic']]):.1f}")

    b = boundaries(grid, plate, W, bk["boundary_threshold_cm_yr"])
    cc = (b["pa"] < n_c) & (b["pb"] < n_c)
    co = (b["pa"] < n_c) != (b["pb"] < n_c)
    oo = (b["pa"] >= n_c) & (b["pb"] >= n_c)
    names = {1: "divergent", -1: "convergent", 0: "transform"}
    for lab, m in (("continent-continent (rifts)", cc), ("continent-ocean (margins)", co), ("ocean-ocean", oo)):
        tot = m.sum()
        print(f"  {lab}: " + ", ".join(f"{names[k]} {np.mean(b['kind'][m] == k):.2f}" for k in (1, -1, 0))
              + f"; mean normal speed {b['vn'][m].mean():+.1f} cm/yr ({tot} edges)")

    # exterior ocean floor: distance to divergent ocean-ocean boundaries / half rate
    ridge = oo & (b["kind"] == 1)
    rc = b["cell"][ridge]
    half = 0.5 * b["vn"][ridge]
    tree = cKDTree(xyz.reshape(-1, 3)[rc])
    ocean = ~cont
    d, idx = tree.query(xyz.reshape(-1, 3)[ocean.reshape(-1)])
    dist_km = 2 * np.arcsin(np.clip(d / 2, 0, 1)) * cs.R_KM
    age = dist_km / (half[idx] / CM_PER_YR)           # km / (km/Myr)
    age = np.minimum(age, bk["max_ocean_age_myr"])
    floor = np.full(cont.shape, np.nan)
    floor[ocean] = age
    oa = grid.area[ocean]
    print(f"exterior ocean floor at breakup: mean age {np.average(age, weights=oa):.0f} Myr, median "
          f"{np.median(age):.0f}, share at the cap {oa[age >= bk['max_ocean_age_myr']].sum() / oa.sum():.2f}")

    np.savez_compressed(os.path.join(work, "breakup.npz"), plate=plate, n_continental=n_c, euler=W,
                        ocean_age_myr=floor, block_plate=block_plate, centre=centre,
                        b_cell=b["cell"], b_other=b["other"], b_kind=b["kind"], b_vn=b["vn"])
    preview(grid, plate, n_c, W, b, floor, os.path.join(work, "E02_breakup.png"))


def preview(grid, plate, n_c, W, b, floor, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Wd, H = 1440, 720
    lon = np.radians(np.linspace(-180, 180, Wd, endpoint=False) + 0.125)
    lat = np.radians(np.linspace(90, -90, H, endpoint=False) - 0.125)
    LO, LA = np.meshgrid(lon, lat)
    f, j, i = cs.latlon_to_cell(LA, LO, grid.n)
    pl = plate[f, j, i]
    rng = np.random.default_rng(5)
    pal = rng.uniform(0.45, 0.95, (plate.max() + 1, 3))
    img = pal[pl].copy()
    oc = pl >= n_c
    age = floor[f, j, i]
    shade = 1 - 0.55 * np.clip(age / 200, 0, 1)
    img[oc] = np.array([0.25, 0.45, 0.8]) * shade[oc][:, None] + 0.15 * pal[pl[oc]]
    fig, ax = plt.subplots(figsize=(14.4, 7.8), dpi=100)
    ax.imshow(img, extent=(-180, 180, -90, 90))
    xyz = grid.xyz.reshape(-1, 3)
    col = {1: "#ffdd00", -1: "#d01010", 0: "#20c020"}
    for k in (1, -1, 0):
        m = b["kind"] == k
        la, lo = np.degrees(cs.xyz_to_latlon(*xyz[b["cell"][m]].T))
        ax.scatter(lo, la, s=0.6, c=col[k], lw=0)
    # velocity arrows on a coarse lat/lon lattice
    gl = np.radians(np.arange(-75, 76, 15))
    go = np.radians(np.arange(-172.5, 180, 15))
    GO, GL = np.meshgrid(go, gl)
    ff, jj, ii = cs.latlon_to_cell(GL, GO, grid.n)
    p = plate[ff, jj, ii]
    r = np.stack([np.cos(GL) * np.cos(GO), np.cos(GL) * np.sin(GO), np.sin(GL)], -1)
    v = np.cross(W[p], r) * cs.R_KM * CM_PER_YR
    east = np.stack([-np.sin(GO), np.cos(GO), np.zeros_like(GO)], -1)
    north = np.stack([-np.sin(GL) * np.cos(GO), -np.sin(GL) * np.sin(GO), np.cos(GL)], -1)
    ax.quiver(np.degrees(GO), np.degrees(GL), (v * east).sum(-1), (v * north).sum(-1), scale=150, width=0.0018)
    ax.set_title("E02 plates at 100 Ma: velocities (arrows), boundaries divergent (yellow), convergent (red), "
                 "transform (green); ocean floor darker = older")
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 30))
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/fw")
