"""E03: 100 Myr of drift after the breakup (fictional_world/tectonics/METHOD.md, step 3).

Semi-Lagrangian rigid-plate advection on the cubed sphere. Each step, every cell looks back along the rotation of each
plate for the cell it came from; a plate claims the cell if that source cell belonged to it.
  * no claim: a gap opened between diverging plates -> new ocean floor born at this time; the cell joins the plate of
    its neighbours (so a ridge is shared between the two plates);
  * one claim: the crust moves with its plate;
  * several claims: plates converge.
      - continent over ocean: the ocean floor is subducted; a continental arc grows behind the trench;
      - ocean over ocean: the older (denser) floor is subducted under the younger; an island arc grows;
      - continent against continent: collision; an orogen; the two plates weld after a sustained collision
        (rotation of the merged plate = area-weighted mean).
Continental cells next to new ocean floor become rifted margins (first time only).

Every continental cell carries its terrane and the age and kind of its last tectonothermal event (codes in
data/tectonics.json, drift.event_codes); every oceanic cell its birth age.

Usage: python3 drift.py WORK_DIR   (reads supercontinent.npz and breakup.npz)
"""

import os
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from supercontinent import CFG, Grid, cs, noise  # noqa: E402
from breakup import boundaries, euler_poles, velocity  # noqa: E402

BASEMENT, RIFT, ARC, COLLISION, ISLAND_ARC, SUTURE, INNER_SUTURE, HOTSPOT, LIP = range(9)


def rotate(p, w, dt):
    """Rotate unit vectors p (N, 3) by the rotation vector w (rad/Myr) over dt Myr (Rodrigues)."""
    ang = np.linalg.norm(w) * dt
    if ang == 0:
        return p
    k = w / np.linalg.norm(w)
    c, s = np.cos(ang), np.sin(ang)
    return p * c + np.cross(k, p) * s + np.outer(p @ k, k) * (1 - c)


def cell_of(grid, p):
    lat, lon = cs.xyz_to_latlon(p[:, 0], p[:, 1], p[:, 2])
    f, j, i = cs.latlon_to_cell(lat, lon, grid.n)
    return (f * grid.n + j) * grid.n + i


def band(grid, src_cells, lo_km, hi_km, among, with_dist=False):
    """Flat indices of cells in `among` whose distance to the nearest of src_cells lies in [lo, hi] (and the
    distances in km if with_dist)."""
    if len(src_cells) == 0:
        return (np.zeros(0, int), np.zeros(0)) if with_dist else np.zeros(0, int)
    xyz = grid.xyz.reshape(-1, 3)
    cand = np.flatnonzero(among)
    tree = cKDTree(xyz[src_cells])
    d, _ = tree.query(xyz[cand], distance_upper_bound=hi_km / cs.R_KM * 1.001)
    ok = np.isfinite(d)
    km = 2 * np.arcsin(np.clip(d[ok] / 2, 0, 1)) * cs.R_KM
    keep = km >= lo_km
    return (cand[ok][keep], km[keep]) if with_dist else cand[ok][keep]


def make_plumes(rng, t0):
    """Mantle plumes in the simulation (mantle) frame: position, slow drift, primary or not, start and end (Ma), LIP
    area for primary plumes (their heads)."""
    pc = CFG["plumes"]
    pts = []
    cosmin = np.cos(np.radians(pc["min_spacing_deg"]))
    while len(pts) < pc["n_primary"] + pc["n_secondary"]:
        v = rng.standard_normal(3)
        v /= np.linalg.norm(v)
        if all(v @ q < cosmin for q in pts):
            pts.append(v)
    plumes = []
    for k, v in enumerate(pts):
        primary = k < pc["n_primary"]
        if primary:
            start = t0 - rng.uniform(0, pc["breakup_window_myr"]) if k < pc["primary_at_breakup"] \
                else rng.uniform(0, t0 - pc["breakup_window_myr"])
            end = -1.0
            lo, hi = np.log(pc["lip_area_km2"])
            lip = float(np.exp(rng.uniform(lo, hi)))
        else:
            start = rng.uniform(0, t0)
            end = start - rng.uniform(*pc["secondary_duration_myr"])
            lip = 0.0
        plumes.append(dict(xyz=v, primary=primary, start=float(start), end=float(end), lip_km2=lip))
    # slow drift of the plumes in the mantle: a random direction and a speed up to max_drift_cm_yr
    for pl in plumes:
        d = rng.standard_normal(3)
        d -= (d @ pl["xyz"]) * pl["xyz"]
        d /= np.linalg.norm(d)
        spd = rng.uniform(0, pc["max_drift_cm_yr"]["primary" if pl["primary"] else "secondary"])
        pl["w"] = np.cross(pl["xyz"], d) * (spd / 0.1 / cs.R_KM)       # rad/Myr
    return plumes


def _pair_sums(pa, pb, w):
    out = {}
    for p, q, x in zip(np.minimum(pa, pb), np.maximum(pa, pb), w):
        out[(int(p), int(q))] = out.get((int(p), int(q)), 0.0) + x
    return out


def main(work):
    s = dict(np.load(os.path.join(work, "supercontinent.npz")))
    b = dict(np.load(os.path.join(work, "breakup.npz")))
    dcfg = CFG["drift"]
    t0 = CFG["breakup_age_ma"]
    dt = dcfg["step_myr"]
    grid = Grid(CFG["sim_n"])
    xyz = grid.xyz.reshape(-1, 3)
    area = grid.area.reshape(-1)
    nc = xyz.shape[0]

    # ------------------------------------------------------------ state at the breakup
    plate = b["plate"].reshape(-1).copy()
    W = b["euler"].copy()
    n_p = len(W)
    cont = s["continent"].reshape(-1).copy()
    terrane = s["terrane"].reshape(-1).copy()
    birth = np.where(cont, np.nan, t0 + b["ocean_age_myr"].reshape(-1))      # Ma
    sage = s["suture_age_ma"].reshape(-1)
    between = s["suture_between_blocks"].reshape(-1)
    tage = s["terrane_age_ma"]
    ev_age = np.where(cont, tage[np.maximum(terrane, 0)], np.nan)
    ev_kind = np.where(cont, BASEMENT, -1)
    insut = cont & np.isfinite(sage)
    ev_age[insut] = sage[insut]
    ev_kind[insut] = np.where(between[insut], SUTURE, INNER_SUTURE)
    # crustal thickness (km): cratons, other basement and old sutures of the supercontinent, ocean crust
    tk = dcfg["thickness_km"]
    crat0 = np.where(terrane >= 0, s["terrane_craton"][np.maximum(terrane, 0)], False)
    th = np.where(cont, np.where(crat0, tk["craton"], tk["basement"]), tk["ocean"]).astype(float)
    th[insut] = tk["old_orogen"]
    vol0 = (area * th)[cont].sum()
    vol_eroded = 0.0
    plumes = make_plumes(np.random.default_rng(CFG["seed"] + 6), t0)
    plc = CFG["plumes"]
    # volcano centres along the tracks: position, plate, birth (Ma), primary plume, born on continental crust
    vol = dict(xyz=np.zeros((0, 3)), plate=np.zeros(0, int), birth=np.zeros(0), primary=np.zeros(0, bool),
               cont=np.zeros(0, bool))
    vol_dist = np.zeros(len(plumes))
    vol_lost = [0, 0]                      # subducted, overridden by continent
    track_age = np.full(nc, np.nan)       # Ma when the crust passed over a plume
    lip = np.zeros(nc, bool)
    short_debt = np.zeros(n_p)             # km2 of continental area each plate owes to post-collisional shortening
    shortened = [0.0]
    accreted = [0.0]
    coll_until = np.full(nc, -np.inf)     # Ma until which a welded collision zone keeps converging
    arc_time = np.zeros(nc)               # Myr a cell of ocean floor has spent in an island-arc band
    fixed_area = [0.0, 0.0]               # area added / removed by the conservative correction
    rifted = np.zeros(nc, bool)
    island = np.zeros(nc, bool)          # island-arc crust on the ocean floor
    collide_count = {}
    merges = []
    cont_area0 = area[cont].sum()

    snaps = {}
    log = []

    def snapshot(t):
        snaps[t] = dict(lip=lip.copy(), plate=plate.copy(), cont=cont.copy(), th=th.copy(), birth=birth.copy(), ev_kind=ev_kind.copy(),
                        ev_age=ev_age.copy(), island=island.copy())

    snapshot(t0)
    tic = time.time()
    rng_pl = np.random.default_rng(CFG["seed"] + 7)
    lips = []
    welds = []
    crat = np.where(terrane >= 0, s["terrane_craton"][np.maximum(terrane, 0)], False)
    rng = np.random.default_rng(CFG["seed"] + 3)
    reorgs = []
    for step in range(int(t0 / dt)):
        t = t0 - (step + 1) * dt          # time after this step, Ma
        if (t + dt) in dcfg["reorganisation_ma"]:
            # plate reorganisation: new Euler poles from the present configuration
            crat = cont & np.where(terrane >= 0, s["terrane_craton"][np.maximum(terrane, 0)], False)
            Wn, _ = euler_poles(grid, plate.reshape(grid.lat.shape), cont.reshape(grid.lat.shape),
                                crat.reshape(grid.lat.shape), CFG["breakup"], rng)
            live = np.unique(plate)
            W[live] = Wn[live]
            reorgs.append(t + dt)
        live = np.unique(plate)
        parea = np.bincount(plate, weights=area, minlength=n_p)
        a_before = np.bincount(plate[cont], weights=area[cont], minlength=n_p)
        lost_by_plate = np.zeros(n_p)
        acc_cells, acc_area, acc_th = [], [], []

        claim = np.zeros((n_p, nc), bool)
        src = np.zeros((n_p, nc), int)
        for p in live:
            sp = cell_of(grid, rotate(xyz, W[p], -dt))
            src[p] = sp
            claim[p] = plate[sp] == p
        n_claim = claim.sum(0)

        # winner among the claims
        score = np.full((n_p, nc), -np.inf)
        for p in live:
            sp = src[p]
            c = cont[sp]
            # continental crust beats ocean; between continents the larger plate overrides; between oceans the
            # younger floor (larger birth age in Ma = older, so use -birth)
            sc = np.where(c, 1e6 + parea[p] / 1e6, -birth[sp])
            score[p] = np.where(claim[p], sc, -np.inf)
        win = score.argmax(0)
        has = n_claim > 0
        s_win = src[win, np.arange(nc)]

        new_plate = np.where(has, win, -1)
        new_cont = np.where(has, cont[s_win], False)
        new_terr = np.where(has, terrane[s_win], -1)
        new_birth = np.where(has, birth[s_win], np.nan)
        new_ev_age = np.where(has, ev_age[s_win], np.nan)
        new_ev_kind = np.where(has, ev_kind[s_win], -1)
        new_rifted = np.where(has, rifted[s_win], False)
        new_island = np.where(has, island[s_win], False)
        new_th = np.where(has, th[s_win], tk["ocean"])
        new_track = np.where(has, track_age[s_win], np.nan)
        new_lip = np.where(has, lip[s_win], False)
        new_arc_time = np.where(has, arc_time[s_win], 0.0)
        new_coll = np.where(has, coll_until[s_win], -np.inf)

        # convergence: what lost
        multi = np.flatnonzero(n_claim > 1)
        lost_cont = np.zeros(nc, bool)
        lost_ocean = np.zeros(nc, bool)
        lost_vol = np.zeros(nc)           # continental crust (km3) overridden at each cell
        for p in live:
            m = multi[(claim[p, multi]) & (win[multi] != p)]
            c = cont[src[p, m]]
            lost_cont[m[c]] = True
            lost_ocean[m[~c]] = True
            np.add.at(lost_vol, m[c], area[src[p, m[c]]] * th[src[p, m[c]]])
            lost_by_plate[p] += area[src[p, m[c]]].sum()
            # buoyant ocean crust (island arcs, plateaus) meeting a continent is accreted, not subducted
            mo = m[~c]
            so = src[p, mo]
            thick = (th[so] >= dcfg["accretion_min_thickness_km"]) & cont[s_win[mo]]
            acc_cells.append(mo[thick])
            acc_area.append(area[so[thick]])
            acc_th.append(th[so[thick]])
        cont_lost_area = area[lost_cont & new_cont].sum()

        # gaps: new ocean floor
        gap = ~has
        new_cont[gap] = False
        new_birth[gap] = t
        new_island[gap] = False
        new_th[gap] = tk["ocean"]
        new_track[gap] = np.nan
        new_lip[gap] = False
        new_arc_time[gap] = 0.0
        new_coll[gap] = -np.inf
        # assign gap cells to the plate of their neighbours (majority), growing inwards
        todo = np.flatnonzero(gap)
        for _ in range(50):
            if todo.size == 0:
                break
            nb = new_plate[grid.nb4[todo]]
            ok = (nb >= 0).any(1)
            for idx in np.flatnonzero(ok):
                v = nb[idx][nb[idx] >= 0]
                new_plate[todo[idx]] = np.bincount(v).argmax()
            todo = todo[~ok]

        plate, cont, terrane, birth = new_plate, new_cont, new_terr, new_birth
        ev_age, ev_kind, rifted, island, th = new_ev_age, new_ev_kind, new_rifted, new_island, new_th
        track_age, lip, arc_time, coll_until = new_track, new_lip, new_arc_time, new_coll

        # nearest-cell resampling does not conserve area: a rotated grid skips some source cells and repeats others,
        # and at plate edges the losses do not cancel. Each plate's continental area is restored to its value before
        # the step minus the crust overridden in collisions (a conservative correction, as in mass-conserving
        # semi-Lagrangian schemes): missing area is added at the plate's own coast, filling the cells with the most
        # continental neighbours first; surplus is removed from the most exposed coastal cells.
        a_after = np.bincount(plate[cont], weights=area[cont], minlength=n_p)
        deficit = a_before - lost_by_plate - a_after - short_debt
        short_debt[:] = 0
        nbp = plate[grid.nb4]
        nbc = cont[grid.nb4]
        for p in np.flatnonzero(np.abs(deficit) > 0.5 * area.mean()):
            same = (nbp == p) & nbc
            k = same.sum(1)
            if deficit[p] > 0:
                cand = np.flatnonzero(~cont & (plate == p) & (k > 0))
                if cand.size == 0:
                    continue
                cand = cand[np.argsort(-k[cand], kind="stable")]
                take = cand[:np.searchsorted(np.cumsum(area[cand]), deficit[p]) + 1]
                src_nb = grid.nb4[take, np.argmax(same[take], axis=1)]
                cont[take] = True
                terrane[take] = terrane[src_nb]
                ev_age[take], ev_kind[take] = ev_age[src_nb], ev_kind[src_nb]
                th[take] = th[src_nb]
                birth[take] = np.nan
                island[take] = False
                rifted[take] = rifted[src_nb]
                fixed_area[0] += area[take].sum()
            else:
                cand = np.flatnonzero(cont & (plate == p) & (k < 4))
                # spare active arcs, so that a correction does not make the margin retreat and the arc migrate
                recent_arc = (ev_kind[cand] == ARC) & (ev_age[cand] <= t + dcfg["spare_arc_myr"])
                cand = np.concatenate([cand[~recent_arc][np.argsort(k[cand[~recent_arc]], kind="stable")],
                                       cand[recent_arc][np.argsort(k[cand[recent_arc]], kind="stable")]])
                take = cand[:np.searchsorted(np.cumsum(area[cand]), -deficit[p]) + 1]
                cont[take] = False
                terrane[take] = -1
                ev_age[take], ev_kind[take] = np.nan, -1
                th[take] = tk["ocean"]
                birth[take] = t
                fixed_area[1] += area[take].sum()

        # accretion: the accreted crust is added in front of the continent, at the ocean cells nearest to where it
        # arrived that touch the continent of the winning plate
        if acc_cells and sum(len(x) for x in acc_cells):
            ac = np.concatenate(acc_cells)
            aa = np.concatenate(acc_area)
            at = np.concatenate(acc_th)
            for w in np.unique(plate[ac]):
                sel = plate[ac] == w
                need = aa[sel].sum()
                touch = ((plate[grid.nb4] == w) & cont[grid.nb4]).any(1)
                cand = np.flatnonzero(~cont & touch)
                if cand.size == 0:
                    continue
                tree = cKDTree(xyz[ac[sel]])
                dd, _ = tree.query(xyz[cand])
                ok = dd * cs.R_KM <= dcfg["accretion_reach_km"]
                cand = cand[ok][np.argsort(dd[ok], kind="stable")]
                take = cand[:np.searchsorted(np.cumsum(area[cand]), need) + 1] if cand.size else cand
                if take.size == 0:
                    continue
                nb = grid.nb4[take]
                srcn = nb[np.arange(take.size), np.argmax((plate[nb] == w) & cont[nb], axis=1)]
                cont[take] = True
                plate[take] = w
                terrane[take] = terrane[srcn]
                th[take] = np.average(at[sel], weights=aa[sel])
                ev_age[take], ev_kind[take] = t, ARC
                birth[take] = np.nan
                island[take] = False
                accreted[0] += area[take].sum()

        # mantle plumes: a plume head arriving makes a large igneous province; an active plume leaves a track
        # volcano centres ride on their plates; those whose cell has passed to another plate are gone (subducted or
        # overridden)
        if len(vol["birth"]):
            for p in np.unique(vol["plate"]):
                m = vol["plate"] == p
                vol["xyz"][m] = rotate(vol["xyz"][m], W[p], dt)
            under = cell_of(grid, vol["xyz"])
            gone = plate[under] != vol["plate"]
            vol_lost[0] += int((gone & ~cont[under]).sum())
            vol_lost[1] += int((gone & cont[under]).sum())
            for k in vol:
                vol[k] = vol[k][~gone]
        for pl in plumes:
            pl["xyz"] = rotate(pl["xyz"][None], pl["w"], dt)[0]
        for k, pl in enumerate(plumes):
            if not (pl["end"] < t <= pl["start"]):
                continue
            if pl["primary"] and pl["start"] - dt < t <= pl["start"]:
                r = np.sqrt(pl["lip_km2"] / np.pi)
                # an irregular outline: the radius varies by +-lip_outline_noise with a smooth noise field
                wob = noise(grid, rng_pl, 25, 2, 2).reshape(-1)
                dist = np.arccos(np.clip(xyz @ pl["xyz"], -1, 1)) * cs.R_KM
                reff = r * (1 + plc["lip_outline_noise"] * wob)
                cells = np.flatnonzero(dist <= reff)
                # lavas thin towards the rim of the province
                taper = np.clip((reff[cells] - dist[cells]) / plc["lip_rim_taper_km"], 0, 1)
                oc_ = cells[~cont[cells]]
                co_ = cells[cont[cells]]
                lip[cells] = True
                tpl = rng_pl.uniform(*plc["plateau_thickness_km"])
                tap_o = taper[~cont[cells]]
                th[oc_] = np.maximum(th[oc_], tk["ocean"] + (tpl - tk["ocean"]) * tap_o)
                th[co_] = np.minimum(th[co_] + plc["flood_basalt_added_km"] * taper[cont[cells]], tk["max"])
                ev_age[cells] = t
                ev_kind[cells] = np.where(cont[cells], LIP, ev_kind[cells])
                lips.append((t, round(float(np.degrees(np.arcsin(pl["xyz"][2]))), 1), round(pl["lip_km2"] / 1e6, 2),
                             round(float(cont[cells].mean()), 2)))
            kind = "primary" if pl["primary"] else "secondary"
            # new volcano centres every spacing_km of plate motion over the plume
            c0 = cell_of(grid, pl["xyz"][None])[0]
            p0 = plate[c0]
            vrel = velocity(W[p0], pl["xyz"]) - velocity(pl["w"], pl["xyz"])
            vkm = np.linalg.norm(vrel) * 10 * dt            # km this step
            vol_dist[k] += vkm
            sp = plc["volcano_spacing_km"][kind]
            new = []
            while vol_dist[k] >= sp:
                vol_dist[k] -= sp
                # the earlier of this step's volcanoes have already moved downstream with the plate
                frac = vol_dist[k] / max(vkm, 1e-9)
                new.append(rotate(pl["xyz"][None], W[p0], frac * dt)[0])
            if new:
                nn = len(new)
                vol["xyz"] = np.vstack([vol["xyz"], np.array(new)])
                vol["plate"] = np.concatenate([vol["plate"], np.full(nn, p0)])
                vol["birth"] = np.concatenate([vol["birth"], t + dt * np.linspace(0, 1, nn, endpoint=False)[::-1]])
                vol["primary"] = np.concatenate([vol["primary"], np.full(nn, pl["primary"])])
                vol["cont"] = np.concatenate([vol["cont"], np.full(nn, bool(cont[c0]))])
            cells = np.flatnonzero(np.arccos(np.clip(xyz @ pl["xyz"], -1, 1)) * cs.R_KM <= plc["track_radius_km"][kind])
            oc_ = cells[~cont[cells]]
            co_ = cells[cont[cells]]
            th[oc_] = np.maximum(th[oc_], plc["track_thickness_km"][kind])
            track_age[cells] = t
            ev_age[co_] = t
            ev_kind[co_] = np.where(ev_kind[co_] == LIP, LIP, HOTSPOT)

        # overridden continental crust is not destroyed: it thickens the continental crust within the collision band
        # (volume conserved up to max_km; the excess counts as eroded)
        lc = np.flatnonzero(lost_vol > 0)
        if lc.size:
            tgt, dkm = band(grid, lc, 0, dcfg["collision_band_km"], cont, with_dist=True)
            if tgt.size:
                tree = cKDTree(xyz[tgt])
                _, near = tree.query(xyz[lc])
                add = np.zeros(nc)
                # each lost cell's volume goes to the band cells near it: spread over the band within the band
                # width, approximated by giving it to its nearest band cell and smoothing over the band below
                np.add.at(add, tgt[near], lost_vol[lc])
                band_area = area[tgt].sum()
                tot = add.sum()
                # spread: half to the nearest cells, half uniformly over the band of this step
                dth = 0.5 * add[tgt] / area[tgt] + 0.5 * tot / band_area
                new = th[tgt] + dth
                over = np.maximum(new - tk["max"], 0)
                vol_eroded += (over * area[tgt]).sum()
                th[tgt] = np.minimum(new, tk["max"])

        # events, from the kinematics of the boundaries (not from the advection, whose nearest-cell sampling also
        # duplicates cells along transform boundaries)
        bd = boundaries(grid, plate.reshape(grid.lat.shape), W, CFG["breakup"]["boundary_threshold_cm_yr"])
        conv = bd["kind"] == -1
        ca, cb = bd["cell"][conv], bd["other"][conv]
        pa, pb, vn = bd["pa"][conv], bd["pb"][conv], bd["vn"][conv]
        conta, contb = cont[ca], cont[cb]
        # overriding side: continent over ocean; younger ocean over older
        a_over = np.where(conta != contb, conta, np.nan_to_num(birth[ca], nan=0) <= np.nan_to_num(birth[cb], nan=0))
        cc = conta & contb
        over = np.where(a_over, ca, cb)
        over_plate = np.where(a_over, pa, pb)
        under_cont = np.where(a_over, contb, conta)
        over_cont = cont[over]
        arc_src = ~cc & over_cont & ~under_cont
        isl_src = ~cc & ~over_cont & ~under_cont
        for pc in np.unique(over_plate[arc_src]):
            cells = np.unique(over[arc_src & (over_plate == pc)])
            arc = band(grid, cells, 0, dcfg["arc_band_km"][1], cont & (plate == pc))
            ev_age[arc] = t
            ev_kind[arc] = ARC
            th[arc] = np.minimum(th[arc] + dcfg["arc_growth_km_per_myr"] * dt, tk["arc_max"])
        lo, hi = dcfg["island_arc_band_km"]
        for po in np.unique(over_plate[isl_src]):
            cells = np.unique(over[isl_src & (over_plate == po)])
            ia = band(grid, cells, lo, hi, ~cont & (plate == po))
            # arc crust builds up only where the arc stays: a cell swept briefly by a migrating trench is not an arc
            arc_time[ia] += dt
            ia = ia[arc_time[ia] >= dcfg["island_arc_min_myr"]]
            island[ia] = True
            th[ia] = np.minimum(th[ia] + dcfg["arc_growth_km_per_myr"] * dt, tk["island_arc_max"])
            ev_age[ia] = t
            ev_kind[ia] = ISLAND_ARC
        if cc.any():
            orog = band(grid, np.unique(np.concatenate([ca[cc], cb[cc]])), 0, dcfg["collision_band_km"], cont)
            ev_age[orog] = t
            ev_kind[orog] = COLLISION
        new_ocean = np.flatnonzero(gap)
        if new_ocean.size:
            rm, dkm = band(grid, new_ocean, 0, dcfg["rift_margin_band_km"], cont & ~rifted, with_dist=True)
            rifted[rm] = True
            # stretched margin: thickness tapers from the unstretched crust inland to the thinned crust at the
            # continent-ocean boundary
            f = np.clip(dkm / dcfg["rift_margin_band_km"], 0, 1)
            th[rm] = np.minimum(th[rm], tk["margin_min"] + (th[rm] - tk["margin_min"]) * f)
            ev_age[rm] = t
            ev_kind[rm] = RIFT

        # sustained collisions weld plates: consumed area = edge length x convergence x dt (each edge counted from
        # both sides, hence / 2)
        edge_km = np.sqrt(area[ca[cc]])
        for (p, q), a_km2 in _pair_sums(pa[cc], pb[cc], edge_km * (-vn[cc]) * 10 * dt / 2).items():
            collide_count[(p, q)] = collide_count.get((p, q), 0) + a_km2
            if collide_count[(p, q)] >= dcfg["merge_after_collision_km2"] and (plate == p).any() and \
                    (plate == q).any():
                # the collision zone keeps converging inside the welded plate (rigid plates cannot show it)
                pair = cc & (np.minimum(pa, pb) == p) & (np.maximum(pa, pb) == q)
                zone = band(grid, np.unique(np.concatenate([ca[pair], cb[pair]])), 0, dcfg["collision_band_km"], cont)
                pcol = dcfg["post_collision"]
                coll_until[zone] = t - pcol["duration_myr"]
                welds.append((t, int(zone.size)))
                ap, aq = area[plate == p].sum(), area[plate == q].sum()
                keep, gone = (p, q) if ap >= aq else (q, p)
                plate[plate == gone] = keep
                vol["plate"][vol["plate"] == gone] = keep
                # a collision removes the slab pull between the two plates: the welded plate takes a new motion by
                # the same rule as at reorganisations (as the India-Asia collision reorganised plate motions)
                crat = cont & np.where(terrane >= 0, s["terrane_craton"][np.maximum(terrane, 0)], False)
                Wn, _ = euler_poles(grid, plate.reshape(grid.lat.shape), cont.reshape(grid.lat.shape),
                                    crat.reshape(grid.lat.shape), CFG["breakup"], rng)
                W[keep] = Wn[keep]
                merges.append((t, int(keep), int(gone)))
                collide_count[(p, q)] = -np.inf

        # post-collisional convergence: the zone stays an active orogen, thickens and widens outwards (plateau growth)
        act = cont & (coll_until <= t) & np.isfinite(coll_until)
        if act.any():
            pcol = dcfg["post_collision"]
            every = max(1, int(round(np.sqrt(area.mean()) / pcol["widening_km_per_myr"])))
            if (step + 1) % every == 0:
                ring = np.unique(grid.nb4[act].reshape(-1))
                ring = ring[cont[ring] & ~act[ring]]
                # a new ring keeps the end time of its neighbours
                nbv = np.where(act[grid.nb4[ring]], coll_until[grid.nb4[ring]], np.inf)
                coll_until[ring] = nbv.min(1)
                act = cont & (coll_until <= t) & np.isfinite(coll_until)
            ev_age[act] = t
            ev_kind[act] = COLLISION
            new = np.minimum(th[act] + pcol["thickening_km_per_myr"] * dt, tk["max"])
            # thickening by shortening conserves volume: the plate loses the equivalent area (next step)
            ia = np.flatnonzero(act)
            dv = area[ia] * (new - th[act])
            np.add.at(short_debt, plate[ia], dv / new)
            shortened[0] += (dv / new).sum()
            th[act] = new

        if (step + 1) % 10 == 0:
            ocean = ~cont
            age = birth[ocean] - t
            lab = grid.components(cont.reshape(grid.lat.shape))
            ids, cnt = np.unique(lab[lab >= 0], return_counts=True)
            big = (np.bincount(lab[lab >= 0], weights=area[lab >= 0])[ids] > 1e6).sum()
            log.append((t, area[cont].sum() / cont_area0, np.average(age, weights=area[ocean]), big,
                        len(np.unique(plate))))
            print(f"t = {t:3.0f} Ma: continental area {log[-1][1]:.3f} of initial, mean ocean age {log[-1][2]:5.1f} "
                  f"Myr, landmasses > 1 M km2: {big}, plates {log[-1][4]}  ({time.time() - tic:.0f} s)", flush=True)
        if (step + 1) * dt % dcfg["snapshot_every_myr"] == 0:
            snapshot(t)
        _ = cont_lost_area
    vol1 = (area * th)[cont].sum()
    print(f"continental crust volume {vol1 / vol0:.3f} of initial (eroded at the thickness cap "
          f"{vol_eroded / vol0:.3f}); thickness mean {np.average(th[cont], weights=area[cont]):.1f} km, "
          f"max {th[cont].max():.0f}")

    print(f"volcano centres: {len(vol['birth'])} kept ({vol['primary'].sum()} primary, {vol['cont'].sum()} on continents),"
          f" {vol_lost[0]} subducted, {vol_lost[1]} overridden by continents")
    print(f"accreted ocean crust {accreted[0] / 1e6:.1f} M km2; post-collisional shortening {shortened[0] / 1e6:.1f} M km2")
    print(f"conservative correction: added {fixed_area[0] / 1e6:.1f} M km2, removed {fixed_area[1] / 1e6:.1f} M km2 "
          f"over 100 Myr")
    print("LIPs (Ma, lat in sim frame, M km2, continental share):", lips)
    print("welds (Ma, cells in the collision zone):", welds)
    print("plate merges (Ma, kept, absorbed):", merges, "; reorganisations at", reorgs)
    np.savez_compressed(os.path.join(work, "drift.npz"), plate=plate, euler=W, continent=cont, terrane=terrane,
                        ocean_birth_ma=birth, thickness_km=th, event_age_ma=ev_age, hotspot_track_ma=track_age, lip=lip,
                        plume_xyz=np.array([p["xyz"] for p in plumes]),
                        plume_primary=np.array([p["primary"] for p in plumes]),
                        plume_start_ma=np.array([p["start"] for p in plumes]),
                        plume_end_ma=np.array([p["end"] for p in plumes]), event_kind=ev_kind, island_arc=island,
                        volcano_xyz=vol["xyz"], volcano_birth_ma=vol["birth"], volcano_primary=vol["primary"],
                        volcano_on_continent=vol["cont"],
                        rifted=rifted, log=np.array(log))
    np.savez_compressed(os.path.join(work, "drift_snapshots.npz"),
                        **{f"{k}_{t:.0f}": v for t, d in snaps.items() for k, v in d.items()})
    preview(grid, snaps, os.path.join(work, "E03_drift.png"))
    age_hist(grid, cont, birth, os.path.join(work, "E03_ocean_age.png"))


def raster(grid, field, W=720, H=360):
    lon = np.radians(np.linspace(-180, 180, W, endpoint=False) + 0.25)
    lat = np.radians(np.linspace(90, -90, H, endpoint=False) - 0.25)
    LO, LA = np.meshgrid(lon, lat)
    f, j, i = cs.latlon_to_cell(LA, LO, grid.n)
    return field.reshape(grid.lat.shape)[f, j, i]


def preview(grid, snaps, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    times = [100, 70, 40, 0]
    fig, axs = plt.subplots(2, 2, figsize=(15, 8.2), dpi=100)
    rng = np.random.default_rng(5)
    pal = rng.uniform(0.45, 0.95, (64, 3))
    kcol = {RIFT: (0.95, 0.85, 0.2), ARC: (0.9, 0.2, 0.15), COLLISION: (0.55, 0.1, 0.5)}
    for ax, t in zip(axs.flat, times):
        d = snaps[t]
        pl = raster(grid, d["plate"])
        c = raster(grid, d["cont"])
        age = raster(grid, d["birth"]) - t
        ek = raster(grid, d["ev_kind"])
        ea = raster(grid, d["ev_age"])
        isl = raster(grid, d["island"])
        img = np.zeros(pl.shape + (3,))
        sh = 1 - 0.6 * np.clip(np.nan_to_num(age) / 200, 0, 1)
        img[~c] = np.array([0.3, 0.55, 0.9]) * sh[~c][:, None]
        img[c] = pal[pl[c] % 64]
        for k, col in kcol.items():
            m = c & (ek == k) & (ea - t <= 30)
            img[m] = 0.35 * img[m] + 0.65 * np.array(col)
        img[isl & ~c] = (0.95, 0.5, 0.2)
        ax.imshow(img, extent=(-180, 180, -90, 90))
        ax.set_title(f"{t} Ma")
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-90, 91, 30))
    fig.suptitle("E03 drift: continents by plate; events of the last 30 Myr: rift margin (yellow), arc (red), "
                 "collision (purple), island arc (orange); ocean darker = older")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def age_hist(grid, cont, birth, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a = grid.area.reshape(-1)
    oc = ~cont
    fig, ax = plt.subplots(figsize=(7, 4), dpi=100)
    ax.hist(birth[oc], bins=np.arange(0, 301, 10), weights=a[oc] / a[oc].sum(), color="#3a6fb0")
    ax.set_xlabel("ocean floor age, Myr")
    ax.set_ylabel("share of ocean floor")
    ax.set_title(f"E03 present ocean floor: mean {np.average(birth[oc], weights=a[oc]):.0f} Myr (Earth ~64)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/home/claude/work/fw")
