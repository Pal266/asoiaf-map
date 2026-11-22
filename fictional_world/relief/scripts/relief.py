"""Layer 4 "relief" of the fictional world: bedrock elevation (cell mean) and expected relief inside each cell on the
10 km grid (experiment E01).

Input (in WORK_DIR, written by fictional_world/tectonics/scripts/output.py):
  tectonics.npz   crust, province, event_age_ma, isostatic_elevation_m (Airy on crustal thickness, swells, sea level)
  boundaries.npz  present plate boundaries (trenches)
  volcanoes.npz   hotspot volcano centres with their ages, volcanic fronts of today's subduction zones
and data/relief.json (in the repository).

Method:
  * continents: the isostatic surface, plus a plain texture by province and ridges and valleys whose amplitude
    grows with elevation (young belts and arcs rough, plateaus smooth on top);
  * continental margins: a shelf to the shelf break over the inner part of the margin, then a slope to the ocean
    floor (where the margin is no wider than max_width_km; wider flooded continent keeps its isostatic depth);
  * ocean floor: the isostatic floor (GDH1 depth, thickened crust, swells) raised by sediment (pelagic + terrigenous
    apron, fitted to GlobSed in the first world), deepened at trenches, with residual depth and abyssal hills;
  * volcanoes: cones on the hotspot volcano centres (summits subsiding with age, as along the Hawaiian chain) and
    on the volcanic fronts of island and continental arcs.

Output (in WORK_DIR): relief.npz [6, 1000, 1000] float32: elevation_m, subgrid_relief_m, sediment_m;
relief_stats.json; E01_relief.png. Deterministic. Bedrock only: no ice, lakes or river incision (later layers).

Usage: python3 relief.py WORK_DIR
"""

import json
import math
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "grid", "scripts"))
import cubesphere as cs  # noqa: E402

P = json.load(open(os.path.join(HERE, "..", "data", "relief.json")))
R = cs.R_KM
N = cs.N


def unit_xyz(n=N):
    out = np.empty((6, n, n, 3))
    for f in range(6):
        lat, lon = cs.cell_latlon(f, n)
        out[f] = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)
    return out


def fractal(xyz, spec, rng):
    """Isotropic unit-variance field from random plane waves, face by face."""
    waves = []
    norm = 0.0
    for lam in spec["wavelengths_km"]:
        amp = (lam / spec["wavelengths_km"][0]) ** spec["hurst"]
        for _ in range(spec["components_per_octave"]):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            waves.append((amp, 2 * np.pi * R / lam * d, rng.uniform(0, 2 * np.pi)))
        norm += amp ** 2 * spec["components_per_octave"] / 2
    out = np.zeros(xyz.shape[:-1], np.float32)
    for f in range(xyz.shape[0]):
        acc = np.zeros(xyz.shape[1:-1])
        for amp, k, ph in waves:
            acc += amp * np.cos(xyz[f] @ k + ph)
        out[f] = acc / math.sqrt(norm)
    return out


def strike_field(field, xyz, halo, smooth_it):
    """Unit tangent vector across strike (along the gradient of a smoothed field) and the gradient magnitude per
    km, from centred differences on the face grid (across face seams through the halo)."""
    f = cs.box_smooth(field, smooth_it, halo)
    F, J, I = halo
    g = f[F, J, I]
    gi = (g[:, 1:-1, 2:] - g[:, 1:-1, :-2]) / 2
    gj = (g[:, 2:, 1:-1] - g[:, :-2, 1:-1]) / 2
    X = xyz[F, J, I]
    ei = (X[:, 1:-1, 2:] - X[:, 1:-1, :-2]) / 2
    ej = (X[:, 2:, 1:-1] - X[:, :-2, 1:-1]) / 2
    # gradient = gi * ei / |ei|^2 + gj * ej / |ej|^2 (the grid is near-orthogonal), per unit of chord length
    gv = gi[..., None] * ei / (ei * ei).sum(-1, keepdims=True) + gj[..., None] * ej / (ej * ej).sum(-1, keepdims=True)
    gv -= (gv * xyz).sum(-1, keepdims=True) * xyz
    mag = np.linalg.norm(gv, axis=-1)
    unit = gv / np.maximum(mag[..., None], 1e-12)
    return unit, mag / R


def fractal_aniso(xyz, spec, rng, across, strength):
    """Unit-variance fractal from random plane waves, each weighted at every point by |cos|^p of the angle between
    its local wave vector and the across-strike direction: wave crests then run along strike. p = exponent x
    strength (0 = isotropic, 1 = full anisotropy)."""
    waves = []
    for lam in spec["wavelengths_km"]:
        amp = (lam / spec["wavelengths_km"][0]) ** spec["hurst"]
        for _ in range(spec["components_per_octave"]):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            waves.append((amp, 2 * np.pi * R / lam * d, rng.uniform(0, 2 * np.pi)))
    out = np.zeros(xyz.shape[:-1], np.float32)
    p = spec["anisotropy_exponent"]
    for f in range(xyz.shape[0]):
        x = xyz[f]
        a = across[f]
        pe = p * strength[f]
        acc = np.zeros(x.shape[:-1])
        wsum = np.zeros(x.shape[:-1])
        for amp, k, ph in waves:
            kt = k - (x @ k)[..., None] * x
            c = np.abs((kt * a).sum(-1)) / np.maximum(np.linalg.norm(kt, axis=-1), 1e-12)
            w = amp * c ** pe
            acc += w * np.cos(x @ k + ph)
            wsum += w * w
        out[f] = acc / np.sqrt(np.maximum(wsum / 2, 1e-12))
    return out


def dist_km(tree_pts, query, upper_km):
    """Great-circle distance (km) from each query point to the nearest tree point, inf beyond upper_km."""
    tree = cKDTree(tree_pts)
    d, i = tree.query(query, distance_upper_bound=2 * math.sin(min(upper_km / R, math.pi) / 2), workers=-1)
    out = np.full(len(query), np.inf)
    ok = np.isfinite(d)
    out[ok] = 2 * np.arcsin(np.clip(d[ok] / 2, 0, 1)) * R
    return out, i


def by_province(prov, table, default=0.0):
    lut = np.full(256, default, np.float32)
    for k, v in table.items():
        lut[int(k)] = v
    return lut[prov]


def main(work):
    t = np.load(os.path.join(work, "tectonics.npz"))
    crust, prov = t["crust"], t["province"]
    cont = crust == 1
    ocean = ~cont
    base = t["isostatic_elevation_m"].astype(np.float64)
    age = np.nan_to_num(t["event_age_ma"], nan=0).astype(np.float64)
    shape = crust.shape
    area = np.broadcast_to(cs.cell_area_km2(N), shape)
    xyz = unit_xyz()
    flat = xyz.reshape(-1, 3)
    halo = cs.halo_index(N, 1)
    rng = np.random.default_rng(P["seed"])

    # ------------------------------------------------------------ ocean floor
    op = P["ocean"]
    cob = cs.boundary_cells(cont, halo)                      # continental cells next to ocean crust
    d_cont, _ = dist_km(flat[cob.reshape(-1)], flat[ocean.reshape(-1)], 3000)
    dc = np.full(shape, np.inf)
    dc[ocean] = d_cont
    sed = np.zeros(shape)
    sed[ocean] = np.minimum(op["pelagic_m_per_myr"] * age[ocean], op["pelagic_max_m"]) + \
        op["apron_m"] * np.exp(-(dc[ocean] / op["apron_km"]) ** 2)
    floor = base + op["sediment_rise_factor"] * sed
    # trenches at today's subduction zones
    b = np.load(os.path.join(work, "boundaries.npz"))
    tr = (b["kind"] == -1) & ~(b["continental"] & b["other_continental"])
    la, lo = np.radians(b["lat"][tr]), np.radians(b["lon"][tr])
    tpts = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    d_tr, _ = dist_km(tpts, flat[ocean.reshape(-1)], 4 * op["trench_km"])
    dtr = np.full(shape, np.inf)
    dtr[ocean] = d_tr
    # the boundary points are 50 km apart: smooth the trench into a continuous furrow
    trench = cs.box_smooth(np.where(ocean, op["trench_m"] * np.exp(-(dtr / op["trench_km"]) ** 2), 0),
                           op["trench_smooth_iterations"], halo)
    floor -= np.where(ocean, trench, 0)
    floor += op["residual"]["rms_m"] * fractal(xyz, op["residual"], rng)
    floor += op["abyssal_hills"]["rms_m"] * fractal(xyz, op["abyssal_hills"], rng)

    # ------------------------------------------------------------ continents
    pn = P["plain_noise"]
    std = cs.box_smooth(by_province(prov, pn["std_m_by_province"]), 5, halo)
    plain = std * fractal(xyz, pn, rng)
    rg = P["ridges"]
    # ridges run along the belts: across-strike direction from the gradient of crustal thickness; full anisotropy
    # where the gradient is steep (belt flanks), fading to isotropic on plateau tops and plains
    across, gmag = strike_field(t["thickness_km"].astype(np.float64), xyz, halo, rg["strike_smooth_iterations"])
    strength = np.clip(gmag / rg["strike_gradient_km_per_km"], 0, 1).astype(np.float32)
    # an isotropic share breaks the ridges into segments (en echelon), as in real belts
    wi = rg["isotropic_share"]
    n_r = np.sqrt(1 - wi) * fractal_aniso(xyz, rg, rng, across, strength)
    iso_spec = dict(rg, components_per_octave=12)
    n_r += np.sqrt(wi) * fractal(xyz, iso_spec, rng)
    rid = 1 - np.abs(n_r)
    rid = (rid - rid.mean()) / rid.std()
    k = cs.box_smooth(by_province(prov, rg["relief_per_metre_by_province"]), 5, halo)
    land_c = base + plain + k * np.maximum(base, 0) * rid

    # continental margins: shelf, then slope down to the ocean floor
    mg = P["margin"]
    sub = cont & (land_c < 0)
    landm = cont & (land_c >= 0)
    q = flat[sub.reshape(-1)]
    d_land, _ = dist_km(flat[landm.reshape(-1)], q, mg["max_width_km"])
    oc_edge = cs.boundary_cells(ocean, halo)                  # ocean cells next to continental crust
    d_oc, i_oc = dist_km(flat[oc_edge.reshape(-1)], q, mg["max_width_km"])
    floor_edge = floor.reshape(-1)[np.flatnonzero(oc_edge.reshape(-1))]
    m = np.isfinite(d_oc) & np.isfinite(d_land)
    fr = np.zeros(len(q))
    fr[m] = d_land[m] / np.maximum(d_land[m] + d_oc[m], 1e-9)
    D = np.where(m, floor_edge[np.minimum(i_oc, len(floor_edge) - 1)], 0)
    sb = mg["shelf_break_m"]
    s = np.clip((fr - mg["shelf_fraction"]) / (1 - mg["shelf_fraction"]), 0, 1)
    smooth = s * s * (3 - 2 * s)
    target = np.where(fr < mg["shelf_fraction"], -sb * fr / mg["shelf_fraction"], -sb + (D + sb) * smooth)
    vals = land_c[sub]
    vals[m] = np.minimum(target[m], -1.0)
    land_c[sub] = vals

    elev = np.where(cont, land_c, floor)

    # ------------------------------------------------------------ volcanoes
    vp = P["volcanoes"]
    v = np.load(os.path.join(work, "volcanoes.npz"))
    cones = []          # (xyz, summit or None, relief above local, radius km)
    hp = vp["hotspot"]
    for lat_, lon_, a, pr in zip(v["hotspot_lat"], v["hotspot_lon"], v["hotspot_birth_ma"], v["hotspot_primary"]):
        kind = "primary" if pr else "secondary"
        cones.append(("hotspot", lat_, lon_, max(float(a), 0.0), kind))
    ap = vp["arc"]
    rng_v = np.random.default_rng(P["seed"] + 1)
    for lat_, lon_, oncont in zip(v["arc_lat"], v["arc_lon"], v["arc_on_continent"]):
        cones.append(("arc", lat_, lon_, 0.0, "continental" if oncont else "island"))
    tree = cKDTree(flat)
    ef = elev.reshape(-1)
    e0 = ef.copy()                 # the surface before any cone: each volcano stands on it, not on its neighbours
    sg_extra = np.zeros(ef.size)
    n_island = {"hotspot": 0, "arc": 0}
    for kind, lat_, lon_, a, sub_kind in cones:
        la, lo = math.radians(lat_), math.radians(lon_)
        c = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])
        _, ic = tree.query(c)
        local = e0[ic]
        oncont = bool(cont.reshape(-1)[ic])
        if kind == "hotspot":
            if oncont:
                if a > hp["continental_young_myr"]:
                    continue
                z, r = local + hp["continental_height_m"], hp["continental_radius_km"]
            else:
                z = hp["summit_m_at_birth"][sub_kind] - hp["subsidence_m_per_sqrt_myr"] * math.sqrt(a)
                z = max(z, local + hp["seamount_height_m"][sub_kind])
                r = hp["radius_km"][sub_kind]
        else:
            if sub_kind == "island":
                z, r = rng_v.uniform(*ap["island_summit_m"]), ap["island_radius_km"]
            else:
                z, r = local + rng_v.uniform(*ap["continental_relief_m"]), ap["continental_radius_km"]
        if z <= local:
            continue
        idx = tree.query_ball_point(c, 2 * math.sin(r / R / 2))
        idx = np.array(idx, int)
        dd = 2 * np.arcsin(np.clip(np.linalg.norm(flat[idx] - c, axis=1) / 2, 0, 1)) * R
        dd = np.maximum(dd, vp["min_cell_distance_km"])
        cone = z - (z - local) * dd / r
        cone = np.maximum(cone, e0[idx])
        ef[idx] = np.maximum(ef[idx], cone)
        sg_extra[idx] = np.maximum(sg_extra[idx], (z - local) / r * 10)
        if z > 0 and local <= 0:
            n_island[kind] += 1
    elev = ef.reshape(shape)

    # ------------------------------------------------------------ sub-grid relief
    sp = P["subgrid"]
    kk = cs.box_smooth(by_province(prov, sp["per_metre_by_province"]), 5, halo)
    subgrid = np.where(cont, sp["land_base_m"] + kk * np.maximum(elev, 0), sp["ocean_m"])
    subgrid = np.maximum(subgrid, sg_extra.reshape(shape))

    np.savez_compressed(os.path.join(work, "relief.npz"), elevation_m=elev.astype(np.float32),
                        subgrid_relief_m=subgrid.astype(np.float32), sediment_m=sed.astype(np.float32))
    stats(work, elev, cont, area, sed, dc, n_island)
    preview(elev, os.path.join(work, "E01_relief.png"))
    zoom(elev, os.path.join(work, "E01_relief_zoom.png"))


def stats(work, elev, cont, area, sed, dc, n_island):
    land = elev > 0
    sea = ~land
    tot = area.sum()
    la = area[land].sum()
    st = {
        "land_share": round(float(la / tot), 3),
        "mean_land_elevation_m": round(float(np.average(elev[land], weights=area[land]))),
        "land_above_1000_2000_3000_m_share_of_land": [round(float(area[elev > h].sum() / la), 3)
                                                      for h in (1000, 2000, 3000)],
        "mean_ocean_depth_m": round(float(-np.average(elev[sea], weights=area[sea]))),
        "shelf_area_above_200m_M_km2": round(float(area[sea & (elev > -200)].sum() / 1e6), 1),
        "ocean_bands_share": {},
        "highest_cell_m": round(float(elev.max())),
        "deepest_cell_m": round(float(elev.min())),
        "volcanic_islands": n_island,
        "land_on_oceanic_crust_km2": round(float(area[land & ~cont].sum())),
        "sediment_margins_m": None,
    }
    bands = [(0, 500), (500, 2500), (2500, 4000), (4000, 5000), (5000, 6000), (6000, 12000)]
    sa = area[sea].sum()
    for lo, hi in bands:
        st["ocean_bands_share"][f"{lo}-{hi}"] = round(float(area[sea & (-elev >= lo) & (-elev < hi)].sum() / sa), 3)
    oc = ~cont
    marg = oc & (dc < 300)
    st["sediment_margins_m"] = round(float(np.average(sed[marg], weights=area[marg])))
    st["sediment_far_ocean_m"] = round(float(np.average(sed[oc & ~marg], weights=area[oc & ~marg])))
    st["sediment_volume_M_km3"] = round(float((sed * area).sum() / 1e3 / 1e6), 1)
    json.dump(st, open(os.path.join(work, "relief_stats.json"), "w"), indent=1)
    print(json.dumps(st, indent=1))


def preview(elev, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    W, H = 2048, 1024
    lon = np.radians(np.linspace(-180, 180, W, endpoint=False) + 180 / W)
    lat = np.radians(np.linspace(90, -90, H, endpoint=False) - 90 / H)
    LO, LA = np.meshgrid(lon, lat)
    f, j, i = cs.latlon_to_cell(LA, LO, N)
    e = elev[f, j, i]
    gy, gx = np.gradient(e)
    shade = np.clip(0.75 + (gx - gy) / 400, 0.35, 1.15)
    img = np.zeros(e.shape + (3,))
    sea = e <= 0
    dep = np.clip(-e / 6000, 0, 1)
    img[sea] = (np.array([0.62, 0.8, 0.95]) * (1 - dep[sea])[:, None] + np.array([0.08, 0.18, 0.42]) * dep[sea][:, None])
    terr = plt.get_cmap("gist_earth")
    img[~sea] = terr(0.38 + 0.6 * np.clip(e[~sea] / 5000, 0, 1) ** 0.6)[:, :3]
    img[~sea] *= shade[~sea][:, None]
    img = np.clip(img, 0, 1)
    fig, ax = plt.subplots(figsize=(20.48, 10.8), dpi=100)
    ax.imshow(img, extent=(-180, 180, -90, 90))
    ax.set_title("E01 relief: bedrock elevation, cell means on the 10 km grid")
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 30))
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def zoom(elev, path):
    """Shaded close-ups of the windows listed in relief.json (preview_windows: lon0, lon1, lat0, lat1)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    wins = P["preview_windows"]
    fig, axs = plt.subplots(1, len(wins), figsize=(8 * len(wins), 8), dpi=100, squeeze=False)
    for ax, (lo0, lo1, la0, la1) in zip(axs[0], wins):
        asp = 1 / np.cos(np.radians((la0 + la1) / 2))
        W = 1000
        H = int(W * (la1 - la0) / (lo1 - lo0) * asp)
        LO, LA = np.meshgrid(np.radians(np.linspace(lo0, lo1, W)), np.radians(np.linspace(la1, la0, H)))
        f, j, i = cs.latlon_to_cell(LA, LO, N)
        z = elev[f, j, i]
        gy, gx = np.gradient(z)
        sh = np.clip(0.75 + (gx - gy) / 500, 0.3, 1.15)
        img = plt.get_cmap("terrain")(0.25 + 0.75 * np.clip(z / 5000, 0, 1))[..., :3] * sh[..., None]
        img[z <= 0] = [0.3, 0.5, 0.8]
        ax.imshow(np.clip(img, 0, 1), extent=(lo0, lo1, la0, la1), aspect=asp)
        ax.set_title(f"lon {lo0}..{lo1}, lat {la0}..{la1}")
    fig.suptitle("E01 relief, close-ups (10 km cells)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/fw")
