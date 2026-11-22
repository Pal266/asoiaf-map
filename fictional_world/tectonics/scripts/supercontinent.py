"""E01: the supercontinent before the breakup (fictional_world/tectonics/METHOD.md, step 1).

All continental crust (41 % of the surface) in one landmass, made of terranes:
  * the terranes are grouped into blocks; each block has a cratonic core (old terranes) and younger accreted terranes,
    welded by old internal sutures;
  * the blocks were welded into the supercontinent by younger sutures, in an order given by a spanning tree of the
    block contacts; the last weld is at assembly_age_ma (~120 Myr before the breakup, as Pangaea).

Deterministic from the seed in data/tectonics.json. Writes an .npz to the work directory (not stored, RULES.md rule 4)
and a preview image.

Usage: python3 supercontinent.py OUT_DIR
"""

import json
import os
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "grid", "scripts"))
import cubesphere as cs  # noqa: E402

CFG = json.load(open(os.path.join(HERE, "..", "data", "tectonics.json")))


# ---------------------------------------------------------------- grid helpers
class Grid:
    def __init__(self, n):
        self.n = n
        lat = np.empty((6, n, n))
        lon = np.empty((6, n, n))
        for f in range(6):
            lat[f], lon[f] = cs.cell_latlon(f, n)
        self.lat, self.lon = lat, lon
        self.xyz = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)
        self.area = np.broadcast_to(cs.cell_area_km2(n), (6, n, n)).copy()
        self.halo = cs.halo_index(n, 1)
        F, J, I = self.halo
        flat = (F * n + J) * n + I
        c = flat[:, 1:-1, 1:-1].reshape(-1)
        self.nb4 = np.stack([flat[:, :-2, 1:-1].reshape(-1), flat[:, 2:, 1:-1].reshape(-1),
                             flat[:, 1:-1, :-2].reshape(-1), flat[:, 1:-1, 2:].reshape(-1)], 1)
        assert (c == np.arange(6 * n * n)).all()
        self.tree = cKDTree(self.xyz.reshape(-1, 3))

    def smooth(self, a, it):
        return cs.box_smooth(a, it, self.halo)

    def components(self, mask):
        """Labels of 4-connected components of mask (flat), -1 outside."""
        m = mask.reshape(-1)
        src = np.repeat(np.arange(m.size), 4)
        dst = self.nb4.reshape(-1)
        keep = m[src] & m[dst]
        g = coo_matrix((np.ones(keep.sum()), (src[keep], dst[keep])), shape=(m.size, m.size))
        _, lab = connected_components(g, directed=False)
        lab = np.where(m, lab, -1)
        return lab


def noise(grid, rng, coarse_n, coarse_smooth, fine_smooth):
    """Smooth random field with unit standard deviation: white noise on a coarse cubed sphere, smoothed, sampled at the
    grid cells and smoothed again."""
    c = rng.standard_normal((6, coarse_n, coarse_n))
    c = cs.box_smooth(c, coarse_smooth, cs.halo_index(coarse_n, 1))
    f, j, i = cs.latlon_to_cell(grid.lat, grid.lon, coarse_n)
    a = grid.smooth(c[f, j, i], fine_smooth)
    return (a - a.mean()) / a.std()


def to_xyz(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def angdist(xyz, p):
    return np.arccos(np.clip(xyz @ p, -1, 1))


def random_point_near(rng, centre, max_deg):
    """Uniform point on the sphere within max_deg of centre."""
    cosmax = np.cos(np.radians(max_deg))
    z = rng.uniform(cosmax, 1)
    phi = rng.uniform(0, 2 * np.pi)
    r = np.sqrt(1 - z * z)
    a = np.array([1.0, 0, 0]) if abs(centre[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(centre, a)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(centre, e1)
    return z * centre + r * (np.cos(phi) * e1 + np.sin(phi) * e2)


# ---------------------------------------------------------------- the landmass
def landmass(grid, rng, sc):
    centre = to_xyz(sc["center_lat"], sc["center_lon"])
    f = np.full(grid.lat.shape, -np.inf)
    for _ in range(sc["n_lobes"]):
        p = random_point_near(rng, centre, sc["lobe_spread_deg"])
        r = np.radians(rng.uniform(*sc["lobe_radius_deg"]))
        f = np.maximum(f, 1 - angdist(grid.xyz, p) / r)
    nz = sc["noise"]
    big = noise(grid, rng, nz["coarse_n"], nz["coarse_smooth"], 2)
    fine = noise(grid, rng, nz["fine_n"], nz["fine_smooth"], 1)
    f = f + nz["amplitude"] * ((1 - nz["fine_weight"]) * big + nz["fine_weight"] * fine)

    target = sc["continental_area_share"] * grid.area.sum()
    # threshold, keep the largest component, fill holes (enclosed sea is continental crust of the interior), and
    # adjust the threshold until the area matches
    lo, hi = f.min(), f.max()
    for _ in range(40):
        th = 0.5 * (lo + hi)
        m = f > th
        lab = grid.components(m)
        big_id = np.bincount(lab[lab >= 0], weights=grid.area.reshape(-1)[lab >= 0]).argmax()
        m = (lab == big_id).reshape(m.shape)
        sea = grid.components(~m)
        sea_id = np.bincount(sea[sea >= 0], weights=grid.area.reshape(-1)[sea >= 0]).argmax()
        m = (sea != sea_id).reshape(m.shape)
        a = grid.area[m].sum()
        if a > target:
            lo = th
        else:
            hi = th
    return m


# ---------------------------------------------------------------- terranes, blocks, sutures
def kmeans_sphere(pts, k, rng, it=50):
    c = pts[rng.choice(len(pts), k, replace=False)]
    for _ in range(it):
        lab = (pts @ c.T).argmax(1)
        for q in range(k):
            if (lab == q).any():
                v = pts[lab == q].sum(0)
                c[q] = v / np.linalg.norm(v)
    return lab, c


def terranes(grid, rng, cont, sc):
    n_t = sc["n_terranes"]
    w = (grid.area * cont).reshape(-1)
    seeds = grid.xyz.reshape(-1, 3)[rng.choice(w.size, n_t, replace=False, p=w / w.sum())]
    # warp cell positions by a smooth noise vector field so that terrane boundaries are not straight
    warp = np.stack([noise(grid, rng, 50, 2, 2) for _ in range(3)], -1) * np.radians(sc["terrane_warp_deg"])
    p = grid.xyz + warp
    p /= np.linalg.norm(p, axis=-1, keepdims=True)
    lab = (p.reshape(-1, 3) @ seeds.T).argmax(1).reshape(cont.shape)
    lab = np.where(cont, lab, -1)
    # a warped Voronoi cell may break into pieces; give stray pieces to the neighbouring terrane they touch most
    for _ in range(5):
        changed = False
        for t in range(n_t):
            comp = grid.components(lab == t)
            ids, cnt = np.unique(comp[comp >= 0], return_counts=True)
            if len(ids) <= 1:
                continue
            for q in ids[cnt < cnt.max()]:
                cells = np.flatnonzero(comp == q)
                nb = lab.reshape(-1)[grid.nb4[cells].reshape(-1)]
                nb = nb[(nb >= 0) & (nb != t)]
                if nb.size:
                    lab.reshape(-1)[cells] = np.bincount(nb).argmax()
                    changed = True
        if not changed:
            break
    # renumber to 0..k-1
    ids = np.unique(lab[lab >= 0])
    remap = -np.ones(n_t, int)
    remap[ids] = np.arange(len(ids))
    lab = np.where(lab >= 0, remap[np.maximum(lab, 0)], -1)
    return lab


def contacts(grid, lab):
    """Pairs of labels that touch, with the number of shared cell edges."""
    a = lab.reshape(-1)
    src = np.repeat(np.arange(a.size), 4)
    dst = grid.nb4.reshape(-1)
    la, lb = a[src], a[dst]
    m = (la >= 0) & (lb >= 0) & (la < lb)
    pairs, cnt = np.unique(np.stack([la[m], lb[m]], 1), axis=0, return_counts=True)
    return {tuple(p): int(c) for p, c in zip(pairs, cnt)}


def assemble(grid, rng, ter, sc):
    n_t = ter.max() + 1
    area = np.bincount(ter[ter >= 0], weights=grid.area[ter >= 0], minlength=n_t)
    centre = np.array([(grid.xyz[ter == t] * grid.area[ter == t][:, None]).sum(0) for t in range(n_t)])
    centre /= np.linalg.norm(centre, axis=1, keepdims=True)
    # blocks: clusters of terranes; the terrane nearest the block's centre, and further ones until craton_share of the
    # block's area, are its craton
    tc = contacts(grid, ter)
    blk, bc = kmeans_sphere(centre, sc["n_blocks"], rng)
    # make every block connected through terrane contacts; a disconnected terrane joins a touching block
    for _ in range(10):
        moved = False
        for b in range(sc["n_blocks"]):
            mem = np.flatnonzero(blk == b)
            if len(mem) <= 1:
                continue
            adj = np.zeros((n_t, n_t), bool)
            for (p, q) in tc:
                adj[p, q] = adj[q, p] = True
            sub = adj[np.ix_(mem, mem)]
            nc, cl = connected_components(coo_matrix(sub), directed=False)
            if nc > 1:
                main = np.bincount(cl, weights=area[mem]).argmax()
                for t in mem[cl != main]:
                    nbb = [blk[q] for q in np.flatnonzero(adj[t]) if blk[q] != b]
                    if nbb:
                        blk[t] = max(set(nbb), key=nbb.count)
                        moved = True
        if not moved:
            break
    ids = np.unique(blk)
    blk = np.searchsorted(ids, blk)
    n_b = len(ids)

    craton = np.zeros(n_t, bool)
    age = np.zeros(n_t)
    for b in range(n_b):
        mem = np.flatnonzero(blk == b)
        c = centre[mem].sum(0)
        c /= np.linalg.norm(c)
        order = mem[np.argsort(-(centre[mem] @ c))]
        acc = 0.0
        for t in order:
            if acc < sc["craton_share"] * area[mem].sum():
                craton[t] = True
                acc += area[t]
        age[mem] = np.where(craton[mem], rng.uniform(*sc["craton_age_ma"], len(mem)),
                            rng.uniform(*sc["terrane_age_ma"], len(mem)))

    # sutures inside a block: the terrane was accreted some time after it formed, before the block was complete
    bcont = {}
    suture_age = {}
    for (p, q), c in tc.items():
        if blk[p] == blk[q]:
            young = min(age[p], age[q])
            suture_age[(p, q)] = max(sc["block_complete_ma"], rng.uniform(0.55, 0.9) * young)
        else:
            key = (min(blk[p], blk[q]), max(blk[p], blk[q]))
            bcont[key] = bcont.get(key, 0) + c
    # assembly of the supercontinent: blocks join along a spanning tree of their contacts (longest contacts first, as
    # the collision closes the ocean between them along the whole length); the welds are dated from
    # block_complete_ma down to assembly_age_ma in a random order
    rows, cols, wts = zip(*[(p, q, 1.0 / c) for (p, q), c in bcont.items()])
    mst = minimum_spanning_tree(coo_matrix((wts, (rows, cols)), shape=(n_b, n_b))).tocoo()
    edges = list(zip(mst.row, mst.col))
    rng.shuffle(edges)
    t_old = sc["assembly_start_ma"]
    ages = np.sort(rng.uniform(sc["assembly_age_ma"], t_old, len(edges)))[::-1]
    ages[-1] = sc["assembly_age_ma"]
    parent = list(range(n_b))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    weld = {}
    for (p, q), a in zip(edges, ages):
        rp, rq = find(p), find(q)
        # every contact between the two growing clusters closes at this collision
        for (u, v) in bcont:
            if (find(u), find(v)) in ((rp, rq), (rq, rp)):
                weld[(u, v)] = a
        parent[rp] = rq
    for (p, q) in tc:
        if blk[p] != blk[q]:
            suture_age[(p, q)] = weld[(min(blk[p], blk[q]), max(blk[p], blk[q]))]
    return dict(block=blk, craton=craton, age=age, suture_age=suture_age, block_weld=weld, n_blocks=n_b)


def suture_belts(grid, ter, asm, width_km):
    """Cells within width/2 of a terrane boundary: the age of that suture and whether it joins two blocks."""
    a = ter.reshape(-1)
    nb = a[grid.nb4]
    edge = (a[:, None] >= 0) & (nb >= 0) & (nb != a[:, None])
    cells, k = np.nonzero(edge)
    p, q = a[cells], nb[cells, k]
    key = np.stack([np.minimum(p, q), np.maximum(p, q)], 1)
    sage = np.array([asm["suture_age"][tuple(x)] for x in key])
    inter = asm["block"][key[:, 0]] != asm["block"][key[:, 1]]
    tree = cKDTree(grid.xyz.reshape(-1, 3)[cells])
    d, idx = tree.query(grid.xyz.reshape(-1, 3), distance_upper_bound=width_km / 2 / cs.R_KM)
    hit = np.isfinite(d) & (a >= 0)
    age = np.where(hit, sage[np.minimum(idx, len(sage) - 1)], np.nan)
    between = np.where(hit, inter[np.minimum(idx, len(sage) - 1)], False)
    return age.reshape(ter.shape), between.reshape(ter.shape)


# ---------------------------------------------------------------- preview
def preview(grid, ter, asm, sage, between, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    W, H = 1440, 720
    lon = np.radians(np.linspace(-180, 180, W, endpoint=False) + 0.125)
    lat = np.radians(np.linspace(90, -90, H, endpoint=False) - 0.125)
    LO, LA = np.meshgrid(lon, lat)
    f, j, i = cs.latlon_to_cell(LA, LO, grid.n)
    t = ter[f, j, i]
    b = np.where(t >= 0, asm["block"][np.maximum(t, 0)], -1)
    cr = np.where(t >= 0, asm["craton"][np.maximum(t, 0)], False)
    rng = np.random.default_rng(3)
    pal = rng.uniform(0.35, 0.9, (asm["n_blocks"], 3))
    img = np.zeros((H, W, 3))
    img[:] = (0.12, 0.2, 0.38)
    land = t >= 0
    img[land] = pal[b[land]] * np.where(cr[land], 0.75, 1.0)[:, None]
    sg = sage[f, j, i]
    bw = between[f, j, i]
    img[np.isfinite(sg) & ~bw] = img[np.isfinite(sg) & ~bw] * 0.6
    img[bw] = (0.85, 0.15, 0.1)
    fig, ax = plt.subplots(figsize=(14.4, 7.6), dpi=100)
    ax.imshow(img, extent=(-180, 180, -90, 90))
    ax.set_title("E01 supercontinent: blocks (colours), cratons (darker), internal sutures (dark), "
                 "assembly sutures (red)")
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 30))
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main(out):
    os.makedirs(out, exist_ok=True)
    sc = CFG["supercontinent"]
    rng = np.random.default_rng(CFG["seed"])
    grid = Grid(CFG["sim_n"])
    cont = landmass(grid, rng, sc)
    ter = terranes(grid, rng, cont, sc)
    asm = assemble(grid, rng, ter, sc)
    sage, between = suture_belts(grid, ter, asm, sc["suture_width_km"])

    total = grid.area.sum()
    print(f"continental crust {grid.area[cont].sum() / 1e6:.1f} M km2 = {grid.area[cont].sum() / total:.3f} of surface")
    n_t = ter.max() + 1
    print(f"terranes {n_t}, blocks {asm['n_blocks']}")
    barea = np.bincount(asm["block"][ter[ter >= 0]], weights=grid.area[ter >= 0])
    for b in range(asm["n_blocks"]):
        mem = np.flatnonzero(asm["block"] == b)
        print(f"  block {b:2d}: {barea[b] / 1e6:5.1f} M km2, {len(mem)} terranes, craton ages "
              f"{', '.join(f'{x:.0f}' for x in asm['age'][mem][asm['craton'][mem]])} Ma")
    welds = sorted(set(round(a) for a in asm["block_weld"].values()), reverse=True)
    print("assembly welds (Ma):", welds)
    crat = asm["craton"][ter[ter >= 0]]
    print(f"cratonic share of continental crust {grid.area[ter >= 0][crat].sum() / grid.area[cont].sum():.2f}")
    print(f"suture belts {grid.area[np.isfinite(sage)].sum() / grid.area[cont].sum():.2f} of continental crust, "
          f"assembly sutures {grid.area[between].sum() / grid.area[cont].sum():.2f}")

    sa = np.array([[p, q, a] for (p, q), a in asm["suture_age"].items()])
    np.savez_compressed(os.path.join(out, "supercontinent.npz"), continent=cont, terrane=ter,
                        terrane_block=asm["block"], terrane_craton=asm["craton"], terrane_age_ma=asm["age"],
                        suture_pairs=sa, suture_age_ma=sage, suture_between_blocks=between)
    preview(grid, ter, asm, sage, between, os.path.join(out, "E01_supercontinent.png"))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/fw")
