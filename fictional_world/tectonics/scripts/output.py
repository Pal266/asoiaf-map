"""E04: the tectonics layer on the 10 km grid (fictional_world/tectonics/METHOD.md, step 4).

From the present state of the drift simulation (drift.npz, 6 x 200 x 200):
  * tectonic provinces in the codes of the first world's history layer (data/tectonics.json, output.province_codes);
  * crust type, age of the last tectonic event (floor age in the oceans), crustal thickness, plate;
  * depth of the ocean floor from its age (GDH1, Stein & Stein 1992);
  * large igneous provinces (oceanic plateaus, flood basalts), hotspot tracks and the swells of active plumes;
  * the spin axis placed where the least land lies near the poles, within a limited polar wander (decision 7): the
    whole planet is rotated from the simulation frame;
  * a first isostatic elevation (Airy on crustal thickness, the continents shifted to Earth's mean height of
    continental crust) and the sea level from the mean basin depth compared with Earth's (water loading 0.69).
    This is a starting surface for the relief layer, not the relief itself.
Fine cells take the coarse cell under a position warped by a smooth fractal field (so that boundaries are not 50 km
staircases); continuous fields are then smoothed.

Usage: python3 output.py WORK_DIR   (reads drift.npz; writes tectonics.npz, boundaries.npz, stats.json, previews)
"""

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from supercontinent import CFG, Grid, cs, noise  # noqa: E402
from breakup import boundaries  # noqa: E402
from drift import ARC, BASEMENT, COLLISION, INNER_SUTURE, RIFT, SUTURE, band  # noqa: E402

R = cs.R_KM


def depth_from_age(t):
    """Basement depth (m) of the ocean floor from its age (Myr), plate model GDH1 (Stein & Stein 1992)."""
    t = np.maximum(t, 0)
    return np.where(t < 20, 2600 + 365 * np.sqrt(t), 5651 - 2473 * np.exp(-0.0278 * t))


def fractal(xyz, wavelengths, hurst, comps, rng):
    """Isotropic unit-variance field from random plane waves (as in the first world's history layer)."""
    f = np.zeros(xyz.shape[:-1])
    norm = 0.0
    for lam in wavelengths:
        amp = (lam / wavelengths[0]) ** hurst
        for _ in range(comps):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            f += amp * np.cos(2 * np.pi * (R / lam) * (xyz @ d) + rng.uniform(0, 2 * np.pi))
        norm += amp ** 2 * comps / 2
    return f / math.sqrt(norm)


def provinces(grid, d, oc, rng):
    """Province codes on the coarse grid."""
    cont = d["continent"]
    kind = d["event_kind"]
    age = d["event_age_ma"]
    th = d["thickness_km"]
    s = np.load(os.path.join(WORK, "supercontinent.npz"))
    ter = d["terrane"]
    crat = cont & np.where(ter >= 0, s["terrane_craton"][np.maximum(ter, 0)], False)
    young = oc["young_ma"]
    pv = np.zeros(cont.shape, np.uint8)
    pv[cont] = 2
    shield = noise(grid, rng, 50, 2, 2).reshape(-1)
    shield = shield > np.quantile(shield[crat], 1 - oc["shield_share_of_cratons"]) if crat.any() else shield > 9
    pv[crat & shield] = 1
    pv[cont & np.isin(kind, [SUTURE, INNER_SUTURE])] = 3
    pv[cont & (kind == BASEMENT) & ~crat & (age < oc["old_orogen_max_ma"])] = 3
    pv[cont & (kind == RIFT)] = 6
    pv[cont & (kind == ARC)] = np.where(age[cont & (kind == ARC)] <= young, 5, 3)
    pv[cont & (kind == COLLISION)] = np.where(age[cont & (kind == COLLISION)] <= young, 4, 3)
    pv[cont & (pv == 4) & (th >= oc["plateau_min_thickness_km"])] = 7
    # foreland and back-arc basins: old basement 100-400 km from young orogens and arcs
    src = np.flatnonzero(cont & np.isin(pv, [4, 5, 7]))
    lo, hi = oc["foreland_band_km"]
    fb = band(grid, src, lo, hi, cont & np.isin(pv, [1, 2]))
    pv[fb] = 8
    pv[~cont & d["island_arc"]] = 5
    pv[~cont & d["lip"]] = 9
    return pv


def choose_frame(grid, d, oc):
    """Rotation matrix (new = M @ old) of the whole planet that places the spin axis, within max_polar_wander_deg of
    the simulation's axis, where the least land lies poleward of polar_latitude_deg (decision 7). Rigid plate motions do not fix the spin axis; true polar wander can
    move it by tens of degrees per 100 Myr. Land is judged on the coarse grid with the same isostasy and sea level as
    the fine output. The rotation is the smallest one that moves the axis (about the axis perpendicular to the old and the new
    spin axis), so the map stays close to the simulation frame."""
    fr = oc["frame"]
    iso = oc["isostasy"]
    cont = d["continent"]
    area = grid.area.reshape(-1)
    k = 1 - iso["rho_crust"] / iso["rho_mantle"]
    e = np.where(cont, 1000 * k * (d["thickness_km"] - iso["reference_thickness_km"]), 0.0)
    e[cont] += iso["mean_continental_elevation_m"] - np.average(e[cont], weights=area[cont])
    depth = depth_from_age(np.nan_to_num(d["ocean_birth_ma"], nan=0))
    sl = oc["sea_level"]
    tt = np.linspace(0, sl["earth_max_age_myr"], 2001)
    pdf = 1 - tt / sl["earth_max_age_myr"]
    sea = sl["water_loading_factor"] * ((depth_from_age(tt) * pdf).sum() / pdf.sum()
                                        - np.average(depth[~cont], weights=area[~cont]))
    land = cont & (e > sea)
    x = grid.xyz.reshape(-1, 3)[land]
    a = area[land]
    # candidate axes: Fibonacci points on the northern hemisphere (an axis and its opposite are the same)
    m = fr["candidates"]
    i = np.arange(m) + 0.5
    z = 1 - i / m
    phi = np.pi * (1 + 5 ** 0.5) * i
    axes = np.stack([np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z], 1)
    s60 = np.sin(np.radians(fr["polar_latitude_deg"]))
    polar = np.zeros(m)
    for c in range(0, m, 500):
        polar[c:c + 500] = ((np.abs(x @ axes[c:c + 500].T) > s60) * a[:, None]).sum(0)
    # only axes within max_polar_wander_deg of the simulation's axis; among the near-minimal ones the closest
    allowed = np.degrees(np.arccos(np.clip(axes[:, 2], -1, 1))) <= fr["max_polar_wander_deg"]
    polar = np.where(allowed, polar, np.inf)
    tol = polar.min() + fr["tie_tolerance"] * a.sum()
    ok = np.flatnonzero(polar <= tol)
    best = ok[np.argmax(np.abs(axes[ok, 2]))]
    ax = axes[best]
    # rotation taking ax to z
    v = np.cross(ax, [0, 0, 1.0])
    sv, cv = np.linalg.norm(v), ax[2]
    if sv < 1e-12:
        M = np.eye(3)
    else:
        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) / sv
        ang = np.arctan2(sv, cv)
        M = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K
    old = ((np.abs(x[:, 2]) > s60) * a).sum() / a.sum()
    lat_ax, lon_ax = np.degrees(cs.xyz_to_latlon(*ax))
    info = {"axis_in_simulation_frame_lat_lon": [round(float(lat_ax), 1), round(float(lon_ax), 1)],
            "polar_wander_deg": round(float(np.degrees(np.arccos(abs(ax[2])))), 1),
            "land_poleward_share_before": round(float(old), 3),
            "land_poleward_share_after": round(float(polar[best] / a.sum()), 3),
            "max_polar_wander_deg": fr["max_polar_wander_deg"]}
    return M, info


def arc_volcanoes(grid, d, bd, ac, rng):
    """Volcanic fronts of today's subduction zones: points arc_distance_km inside the overriding plate from each
    convergent boundary edge (continent over ocean; younger ocean over older), thinned to one per spacing_km."""
    cont = d["continent"]
    birth = np.nan_to_num(d["ocean_birth_ma"], nan=0)
    xyz = grid.xyz.reshape(-1, 3)
    conv = bd["kind"] == -1
    ca, cb, vn = bd["cell"][conv], bd["other"][conv], bd["vn"][conv]
    a_over = np.where(cont[ca] != cont[cb], cont[ca], birth[ca] <= birth[cb])
    sub = ~(cont[ca] & cont[cb])                                # collisions have no arc
    over = np.where(a_over, ca, cb)[sub]
    under = np.where(a_over, cb, ca)[sub]
    vn = vn[sub]
    r = xyz[over]
    nvec = r - xyz[under]
    nvec -= (nvec * r).sum(1, keepdims=True) * r
    nvec /= np.maximum(np.linalg.norm(nvec, axis=1, keepdims=True), 1e-12)
    ang = ac["arc_distance_km"] / cs.R_KM
    p = r * np.cos(ang) + nvec * np.sin(ang)
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    # keep only fronts that land on the overriding plate
    pc = []
    for q in p:
        lat, lon = cs.xyz_to_latlon(*q)
        f, j, i = cs.latlon_to_cell(np.array([lat]), np.array([lon]), grid.n)
        pc.append((f[0] * grid.n + j[0]) * grid.n + i[0])
    pc = np.array(pc, int)
    ok = d["plate"][pc] == d["plate"][over]
    p, vn, pc = p[ok], vn[ok], pc[ok]
    # thin to one volcano per spacing (greedy, in random order)
    order = rng.permutation(len(p))
    keep = []
    tree_pts = []
    lim = 2 * np.sin(ac["spacing_km"] / cs.R_KM / 2)
    for i in order:
        if tree_pts and (np.linalg.norm(np.array(tree_pts) - p[i], axis=1) < lim).any():
            continue
        keep.append(i)
        tree_pts.append(p[i])
    keep = np.array(keep, int)
    return dict(xyz=p[keep] if keep.size else np.zeros((0, 3)), cont=cont[pc[keep]] if keep.size else np.zeros(0, bool),
                vn=-vn[keep] if keep.size else np.zeros(0))


def main(work):
    global WORK
    WORK = work
    oc = CFG["output"]
    d = {k: v.reshape(-1) if v.ndim == 3 else v for k, v in np.load(os.path.join(work, "drift.npz")).items()}
    grid = Grid(CFG["sim_n"])
    rng = np.random.default_rng(CFG["seed"] + 4)
    n0 = grid.n
    shape0 = grid.lat.shape
    cont = d["continent"]
    pv = provinces(grid, d, oc, rng)
    M, frame_info = choose_frame(grid, d, oc)
    print("frame:", frame_info)
    ocean_age = np.where(cont, np.nan, d["ocean_birth_ma"])
    ev_age = np.where(cont, d["event_age_ma"], ocean_age)
    ev_age = np.where(~cont & d["island_arc"], d["event_age_ma"], ev_age)

    # present plate boundaries (coarse)
    bd = boundaries(grid, d["plate"].reshape(shape0), d["euler"], CFG["breakup"]["boundary_threshold_cm_yr"])
    xyz0 = grid.xyz.reshape(-1, 3)
    blat, blon = np.degrees(cs.xyz_to_latlon(*(xyz0[bd["cell"]] @ M.T).T))
    np.savez_compressed(os.path.join(work, "boundaries.npz"), lat=blat.astype(np.float32), lon=blon.astype(np.float32),
                        plate=bd["pa"].astype(np.uint8), other_plate=bd["pb"].astype(np.uint8),
                        kind=bd["kind"].astype(np.int8), normal_cm_yr=bd["vn"].astype(np.float32),
                        tangential_cm_yr=bd["vt"].astype(np.float32),
                        continental=cont[bd["cell"]], other_continental=cont[bd["other"]])

    # volcanoes: hotspot volcano centres from the drift, and the volcanic fronts of today's subduction zones
    vx = d["volcano_xyz"] @ M.T
    vlat, vlon = np.degrees(cs.xyz_to_latlon(*vx.T))
    av = arc_volcanoes(grid, d, bd, oc["arc_volcanoes"], rng)
    alat, alon = np.degrees(cs.xyz_to_latlon(*(av["xyz"] @ M.T).T)) if len(av["xyz"]) else (np.zeros(0), np.zeros(0))
    np.savez_compressed(os.path.join(work, "volcanoes.npz"),
                        hotspot_lat=vlat.astype(np.float32), hotspot_lon=vlon.astype(np.float32),
                        hotspot_birth_ma=d["volcano_birth_ma"].astype(np.float32),
                        hotspot_primary=d["volcano_primary"], hotspot_born_on_continent=d["volcano_on_continent"],
                        arc_lat=alat.astype(np.float32), arc_lon=alon.astype(np.float32),
                        arc_on_continent=av["cont"], arc_convergence_cm_yr=av["vn"].astype(np.float32))
    print(f"volcanoes: {len(vlat)} hotspot centres, {len(alat)} arc volcanoes ({av['cont'].sum()} continental)")

    th0 = d["thickness_km"]
    # ------------------------------------------------------------ fine grid
    n = oc["n"]
    wp = oc["warp"]
    halo = cs.halo_index(n, 1)
    F = {k: np.zeros((6, n, n), t) for k, t in
         (("crust", np.uint8), ("province", np.uint8), ("event_age_ma", np.float32), ("thickness_km", np.float32),
          ("plate", np.uint8), ("ocean_age", np.float32),
          ("hotspot_track_ma", np.float32), ("lip", bool))}
    for f in range(6):
        lat, lon = cs.cell_latlon(f, n)
        x = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)
        w = np.stack([fractal(x, wp["wavelengths_km"], wp["hurst"], wp["components"], np.random.default_rng(
            CFG["seed"] + 50 + k)) for k in range(3)], -1)
        y = x + w * (wp["amplitude_km"] / R)
        y /= np.linalg.norm(y, axis=-1, keepdims=True)
        y = y @ M          # back to the simulation frame (M is orthogonal: old = M.T @ new)
        ff, jj, ii = cs.latlon_to_cell(*cs.xyz_to_latlon(y[..., 0], y[..., 1], y[..., 2]), n0)
        c = (ff * n0 + jj) * n0 + ii
        F["crust"][f] = cont[c]
        F["province"][f] = pv[c]
        F["event_age_ma"][f] = ev_age[c]
        F["thickness_km"][f] = th0[c]
        F["plate"][f] = d["plate"][c]
        F["ocean_age"][f] = np.nan_to_num(ocean_age[c], nan=-1)
        F["hotspot_track_ma"][f] = d["hotspot_track_ma"][c]
        F["lip"][f] = d["lip"][c]
    # smooth continuous fields within each crust type
    for key, it in (("thickness_km", oc["smooth_iterations"]), ("ocean_age", oc["smooth_iterations"])):
        for crust in (0, 1):
            m = F["crust"] == crust
            if key == "ocean_age" and crust == 1:
                continue
            v = np.where(m, F[key], 0).astype(np.float64)
            wsum = cs.box_smooth(m.astype(np.float64), it, halo)
            vs = cs.box_smooth(v, it, halo) / np.maximum(wsum, 1e-9)
            F[key] = np.where(m, vs, F[key]).astype(np.float32)
    oa = np.where(F["crust"] == 0, F["ocean_age"], np.nan)
    F["event_age_ma"] = np.where((F["crust"] == 0) & (F["province"] == 0), oa, F["event_age_ma"]).astype(np.float32)
    depth = depth_from_age(np.nan_to_num(oa, nan=0))

    # ------------------------------------------------------------ first isostatic surface and sea level
    iso = oc["isostasy"]
    k = 1 - iso["rho_crust"] / iso["rho_mantle"]
    e_cont = 1000 * k * (F["thickness_km"] - iso["reference_thickness_km"]) + iso["reference_elevation_m"]
    # island arcs: thickened ocean crust stands higher than the floor around it (Airy, same factor)
    # thickened ocean crust (island arcs, plateaus, hotspot tracks) stands higher than the floor around it (Airy)
    e_ocean = -depth + 1000 * k * np.maximum(F["thickness_km"] - CFG["drift"]["thickness_km"]["ocean"], 0)
    e = np.where(F["crust"] == 1, e_cont, e_ocean).astype(np.float64)
    area = np.broadcast_to(cs.cell_area_km2(n), e.shape)
    # continents: Airy relief about Earth's mean height of continental crust (shelves included)
    cc = F["crust"] == 1
    shift = iso["mean_continental_elevation_m"] - np.average(e[cc], weights=area[cc])
    e[cc] += shift
    # sea level: ocean basins shallower or deeper than Earth's by their age, with water loading (Muller et al. 2008)
    sl = oc["sea_level"]
    tm = sl["earth_max_age_myr"]
    tt = np.linspace(0, tm, 2001)
    pdf = 1 - tt / tm
    d_earth = float((depth_from_age(tt) * pdf).sum() / pdf.sum())
    oc_m = F["crust"] == 0
    d_ours = float(np.average(depth[oc_m], weights=area[oc_m]))
    sea = sl["water_loading_factor"] * (d_earth - d_ours)
    # swells of the active plumes (dynamic topography, moving with the plume, not with the plate)
    sw = CFG["plumes"]["swell"]
    swell = np.zeros(e.shape, np.float32)
    active = d["plume_end_ma"] < 0
    active &= d["plume_start_ma"] >= 0
    pxyz = d["plume_xyz"][active] @ M.T
    prim = d["plume_primary"][active]
    hot = []
    for f in range(6):
        lat, lon = cs.cell_latlon(f, n)
        x = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], -1)
        for v, pr in zip(pxyz, prim):
            dkm = np.arccos(np.clip(x @ v, -1, 1)) * R
            h = sw["height_m"]["primary" if pr else "secondary"]
            swell[f] += h * np.exp(-(dkm / sw["half_width_km"]) ** 2)
    for v, pr in zip(pxyz, prim):
        la, lo = np.degrees(cs.xyz_to_latlon(*v))
        hot.append({"lat": round(float(la), 2), "lon": round(float(lo), 2), "primary": bool(pr)})
    e = e + swell
    elev = (e - sea).astype(np.float32)
    ocean_depth = np.where(F["crust"] == 0, depth + sea, np.nan).astype(np.float32)

    np.savez_compressed(os.path.join(work, "tectonics.npz"), crust=F["crust"], province=F["province"],
                        event_age_ma=F["event_age_ma"], ocean_depth_m=ocean_depth, thickness_km=F["thickness_km"],
                        plate=F["plate"], isostatic_elevation_m=elev, sea_level_m=np.float32(sea),
                        frame_rotation=M, euler_rad_per_myr=(d["euler"] @ M.T),
                        hotspot_track_ma=F["hotspot_track_ma"], lip=F["lip"], hotspot_swell_m=swell)
    json.dump(hot, open(os.path.join(work, "active_hotspots.json"), "w"), indent=1)

    # ------------------------------------------------------------ checks
    tot = area.sum()
    land = elev > 0
    LAT = np.stack([np.degrees(cs.cell_latlon(f, n)[0]) for f in range(6)])
    oc0 = F["crust"] == 0
    st = {"frame": frame_info,
          "active_hotspots": {"primary": int(prim.sum()), "secondary": int((~prim).sum())},
          "oceanic_plateau_share_of_ocean": round(float(area[oc0 & (F["province"] == 9)].sum() / area[oc0].sum()), 4),
          "hotspot_track_share_of_ocean": round(float(area[oc0 & np.isfinite(F["hotspot_track_ma"])].sum()
                                                      / area[oc0].sum()), 4),
          "flood_basalt_share_of_continental_crust": round(float(area[(F["crust"] == 1) & F["lip"]].sum()
                                                                 / area[F["crust"] == 1].sum()), 4),
          "sea_level_m_above_modern_datum": round(float(sea), 1),
          "mean_basement_depth_m": {"earth_reference": round(d_earth), "this_world": round(d_ours)},
          "land_share": round(float(area[land].sum() / tot), 3),
          "continental_crust_share": round(float(area[F["crust"] == 1].sum() / tot), 3),
          "flooded_share_of_continental_crust": round(float(area[(F["crust"] == 1) & ~land].sum()
                                                            / area[F["crust"] == 1].sum()), 3),
          "land_poleward_of_60_share": round(float(area[land & (np.abs(LAT) > 60)].sum() / area[land].sum()), 3),
          "mean_land_elevation_m": round(float(np.average(elev[land], weights=area[land])), 0),
          "mean_ocean_depth_m": round(float(-np.average(elev[~land], weights=area[~land])), 0),
          "province_share_of_continental_crust": {}}
    cc = F["crust"] == 1
    for code, name in oc["province_codes"].items():
        st["province_share_of_continental_crust"][f"{code} {name}"] = round(
            float(area[cc & (F["province"] == int(code))].sum() / area[cc].sum()), 3)
    # landmasses on the coarse grid (fine land sampled at coarse centres)
    land0 = np.zeros(shape0, bool)
    for f in range(6):
        lat, lon = cs.cell_latlon(f, n0)
        ff, jj, ii = cs.latlon_to_cell(lat, lon, n)
        land0[f] = land[ff, jj, ii]
    lab = grid.components(land0)
    a0 = grid.area.reshape(-1)
    ar = np.bincount(lab[lab >= 0], weights=a0[lab >= 0]) / 1e6
    st["landmasses_M_km2"] = [round(float(x), 1) for x in sorted(ar[ar >= 0.5], reverse=True)]
    hist = np.histogram(elev, bins=np.arange(-7000, 7001, 500), weights=area)[0] / tot
    st["hypsometry_500m_bins_from_-7000"] = [round(float(x), 4) for x in hist]
    json.dump(st, open(os.path.join(work, "tectonics_stats.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in st.items() if not k.startswith("hyps")}, indent=1))
    preview(F, elev, os.path.join(work, "E04_tectonics.png"), oc, hot)


def preview(F, elev, path, oc, hot):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    W, H = 1440, 720
    lon = np.radians(np.linspace(-180, 180, W, endpoint=False) + 0.125)
    lat = np.radians(np.linspace(90, -90, H, endpoint=False) - 0.125)
    LO, LA = np.meshgrid(lon, lat)
    n = F["crust"].shape[1]
    f, j, i = cs.latlon_to_cell(LA, LO, n)
    fig, axs = plt.subplots(2, 1, figsize=(13, 13.5), dpi=100)
    e = elev[f, j, i]
    terr = plt.get_cmap("terrain")
    img = np.zeros(e.shape + (3,))
    sea = e <= 0
    img[sea] = np.array([0.1, 0.25, 0.55]) + np.clip(1 + e[sea] / 6000, 0, 1)[:, None] * np.array([0.25, 0.4, 0.4])
    img[~sea] = terr(0.25 + 0.75 * np.clip(e[~sea] / 5000, 0, 1))[:, :3]
    axs[0].imshow(img, extent=(-180, 180, -90, 90))
    axs[0].set_title("E04 first isostatic surface (Airy on crustal thickness; GDH1 ocean floor; sea level from "
                     "basin depth vs Earth). Not the relief layer.")
    cols = ["#26508c", "#7a4f2a", "#c9a46a", "#8c6bb1", "#d7301f", "#fc8d59", "#fee08b", "#a50026", "#d9d98c",
            "#4fd0e0"]
    pv = F["province"][f, j, i]
    axs[1].imshow(pv, cmap=ListedColormap(cols), vmin=-0.5, vmax=9.5, extent=(-180, 180, -90, 90),
                  interpolation="nearest")
    axs[1].set_title("Provinces: " + ", ".join(f"{k} {v.split(' (')[0]}" for k, v in oc["province_codes"].items()),
                     fontsize=8)
    tr = np.isfinite(F["hotspot_track_ma"][f, j, i])
    axs[1].contour(np.linspace(-180, 180, W), np.linspace(90, -90, H), tr.astype(float), levels=[0.5],
                   colors="white", linewidths=0.5)
    for h in hot:
        for ax in axs:
            ax.plot(h["lon"], h["lat"], "o", ms=6 if h["primary"] else 3.5, mfc="#ffea00" if h["primary"] else
                    "white", mec="black", mew=0.6)
    for ax in axs:
        ax.set_xticks(range(-180, 181, 30))
        ax.set_yticks(range(-90, 91, 30))
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


WORK = None

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/fw")
