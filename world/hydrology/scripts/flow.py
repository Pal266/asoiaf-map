"""Water balance and flow routing on any grid (experiment E12).

Water balance (annual, per cell):
- potential evapotranspiration from monthly temperature (Thornthwaite 1948), with the day length of the latitude;
- actual evapotranspiration from the Budyko curve in Fu's form (Fu 1981; Zhang et al. 2004):
  E / P = 1 + phi - (1 + phi^w)^(1/w), phi = PET / P; runoff R = P - E;
- open water (lakes) evaporates at the potential rate.

Routing:
- depressions filled by priority flood with an epsilon slope (Barnes, Lehman & Mulla 2014), from the sea inwards;
- each cell drains to its steepest-descent neighbour (8 neighbours) on the filled surface;
- discharge accumulated downstream in topological order; lake cells add their net balance P - E_open, which can be
  negative, and the flow is never allowed below zero. So a lake whose evaporation exceeds its inflow ends the river:
  a terminal (endorheic) lake; one with surplus spills through its outlet.
- terminal lakes shrink to the lowest cells of their depression whose net evaporation equals the inflow;
- depressions under lake_min_area_km2, or shallower than lake_min_maxdepth_m at their deepest, are taken as filled
  by alluvium; overflowing ones keep only water deeper than lake_overflow_depth_m, the rest is marsh; the floor of
  a terminal depression up to playa_rise_m above its lake is salt flat.
Deterministic.
"""

import heapq

import numpy as np
from numba import njit

SECONDS_PER_YEAR = 365.25 * 86400
MONTH_DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
MID_DOY = np.array([15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349])


def day_length_h(lat_deg, doy):
    decl = 23.44 * np.pi / 180 * np.sin(2 * np.pi * (284 + doy) / 365)
    x = -np.tan(np.radians(lat_deg)) * np.tan(decl)
    return 24 / np.pi * np.arccos(np.clip(x, -1, 1))


def pet_thornthwaite(T_month, lat_deg):
    """Monthly potential evapotranspiration (mm), T_month (12, ...) in C, lat_deg broadcastable to T[0]."""
    Tp = np.maximum(T_month, 0)
    I = ((Tp / 5) ** 1.514).sum(0)
    a = 6.75e-7 * I ** 3 - 7.71e-5 * I ** 2 + 1.792e-2 * I + 0.49239
    Is = np.maximum(I, 1e-6)
    pet = np.zeros(T_month.shape, np.float32)
    for m in range(12):
        L = day_length_h(lat_deg, MID_DOY[m])
        t = Tp[m]
        base = 16 * (10 * t / Is) ** a
        # above 26.5 C Thornthwaite's table (Willmott et al. 1985 fit)
        hot = -415.85 + 32.24 * t - 0.43 * t ** 2
        pet[m] = np.where(t > 26.5, hot, base) * (L / 12) * (MONTH_DAYS[m] / 30)
    return pet


def runoff_budyko(P_year, PET_year, w):
    """Annual runoff (mm) from Fu's form of the Budyko curve."""
    P = np.maximum(P_year, 1e-6)
    phi = PET_year / P
    E = P * (1 + phi - (1 + phi ** w) ** (1 / w))
    return np.maximum(P - np.maximum(E, 0), 0)


def priority_flood(z, sea, nb, eps):
    """Filled surface (flat indices) by priority flood from the sea. z: heights, sea: bool, nb (k, n) neighbours.
    Returns filled heights (float64) and the order in which cells were reached (downstream first)."""
    seeds = np.flatnonzero(sea & (~sea)[nb].any(0))
    return _flood(z.astype(np.float64), sea, nb, seeds, eps)


@njit(cache=True)
def _flood(z, sea, nb, seeds, eps):
    n = z.size
    filled = z.copy()
    done = sea.copy()
    order = np.empty(n, np.int64)
    no = 0
    heap = [(filled[seeds[0]], seeds[0])]
    for i in range(1, seeds.size):
        heap.append((filled[seeds[i]], seeds[i]))
    heapq.heapify(heap)
    while len(heap) > 0:
        h, c = heapq.heappop(heap)
        for k in range(nb.shape[0]):
            q = nb[k, c]
            if done[q]:
                continue
            done[q] = True
            if filled[q] <= h + eps:
                filled[q] = h + eps
            heapq.heappush(heap, (filled[q], np.int64(q)))
            order[no] = q
            no += 1
    return filled, order[:no]


@njit(cache=True)
def receivers(filled, sea, nb, dist):
    """Steepest-descent neighbour on the filled surface; -1 for sea cells."""
    n = filled.size
    rec = np.full(n, -1, np.int64)
    for c in range(n):
        if sea[c]:
            continue
        best = 0.0
        bi = -1
        for k in range(nb.shape[0]):
            q = nb[k, c]
            s = (filled[c] - filled[q]) / dist[k, c]
            if s > best:
                best = s
                bi = q
        rec[c] = bi
    return rec


@njit(cache=True)
def accumulate(order, rec, local):
    """Discharge: local contributions summed downstream, never below zero. order: downstream first."""
    q = local.copy()
    for i in range(order.size - 1, -1, -1):
        c = order[i]
        if q[c] < 0:
            q[c] = 0.0
        r = rec[c]
        if r >= 0:
            q[r] += q[c]
    return q


def depressions(z, filled, land, nb, min_depth):
    """Label connected depressions (filled - z > min_depth) on land; -1 elsewhere."""
    dep = land & (filled - z > min_depth)
    lab = np.full(z.size, -1, np.int64)
    k = 0
    for s in np.flatnonzero(dep):
        if lab[s] >= 0:
            continue
        stack = [s]
        lab[s] = k
        while stack:
            c = stack.pop()
            for q in nb[:, c]:
                if dep[q] and lab[q] < 0:
                    lab[q] = k
                    stack.append(q)
        k += 1
    return lab, k


def route(z, land, P, PET, area_km2, nb, dist, cfg, ice=None, z_dep=None, forced_lake=None):
    """Full routing on flat arrays. z: surface height (m), land: bool, P, PET: annual mm, area_km2, nb/dist: (k, n)
    neighbours and distances (m), ice: bool (ice-covered cells: runoff = P), z_dep: surface for finding lakes if not
    z, forced_lake: bool cells that are lake whatever the depressions say. Returns a dict of flat arrays:
    runoff_mm, discharge_m3s, lake (bool), terminal lake id, receivers, filled."""
    sea = ~land
    filled, order = priority_flood(z, sea, nb, cfg["fill_epsilon_m"])
    rec = receivers(filled, sea, nb, dist)
    R = runoff_budyko(P, PET, cfg["fu_w"])
    if ice is not None:
        R = np.where(ice, P, R)
    to_m3s = area_km2 * 1e6 / 1000 / SECONDS_PER_YEAR          # mm/yr over the cell -> m3/s
    local = np.where(land, R * to_m3s, 0.0)
    # pass 1: depressions treated as land, to find their inflow
    q1 = accumulate(order, rec, local)
    if z_dep is None:
        lab, k = depressions(z, filled, land, nb, cfg["lake_min_depth_m"])
        fd = filled
    else:      # depressions of another surface (Earth test: the unburnt one)
        fd, _ = priority_flood(z_dep, sea, nb, cfg["fill_epsilon_m"])
        lab, k = depressions(z_dep, fd, land, nb, cfg["lake_min_depth_m"])
    lake = np.zeros(z.size, bool)
    net_evap = np.where(land, np.maximum(PET - P, 0) * to_m3s, 0.0)    # open water's loss, m3/s per cell
    gain_open = (P - PET) * to_m3s
    terminal = np.zeros(z.size, bool)
    marsh = np.zeros(z.size, bool)        # sediment-filled part of an overflowing depression: floodplain, marsh
    playa = np.zeros(z.size, bool)        # dry floor of a terminal depression: salt flat
    lakes_info = []
    for d in range(k):
        cells = np.flatnonzero(lab == d)
        zd_ = (z if z_dep is None else z_dep)
        fl_ = (filled if z_dep is None else fd)
        if area_km2[cells].sum() < cfg["lake_min_area_km2"] or (fl_[cells] - zd_[cells]).max() < cfg["lake_min_maxdepth_m"]:
            continue                   # small or shallow: filled by alluvium
        inside = np.zeros(z.size, bool)
        inside[cells] = True
        out = cells[~inside[rec[cells]]]                       # cells that drain out of the depression
        inflow = float(q1[out].sum())
        # the lake surface: depression cells sorted by height; full if the inflow can feed the net evaporation of
        # the whole depression, else the lowest cells until their net evaporation equals the inflow
        srt = cells[np.argsort((z if z_dep is None else z_dep)[cells])]
        cum = np.cumsum(net_evap[srt])
        if cum[-1] <= inflow:
            # an overflowing depression: over geological time its shallow margins fill with sediment and its outlet
            # incises, so only water deeper than lake_overflow_depth_m remains lake (tectonic and glacial basins)
            zd = (z if z_dep is None else z_dep)
            fl = (filled if z_dep is None else fd)
            wet = cells[fl[cells] - zd[cells] > cfg["lake_overflow_depth_m"]]
            marsh[cells] = True
            term = False
        else:
            nwet = int(np.searchsorted(cum, inflow)) + 1
            wet = srt[:nwet]
            # salt flat: the floor around the terminal lake, up to playa_rise_m above its surface (flooded in wet
            # years); higher parts of a closed basin are alluvial plains and fans
            zs = (z if z_dep is None else z_dep)
            level = zs[wet].max()
            playa[cells[zs[cells] <= level + cfg["playa_rise_m"]]] = True
            term = True
        lake[wet] = True
        terminal[wet] = term
        lakes_info.append({"cells": int(len(wet)), "area_km2": round(float(area_km2[wet].sum())),
                           "terminal": term, "inflow_m3s": round(inflow, 1)})
    if forced_lake is not None:                # lakes given from outside (the canon lakes of the map)
        lake |= forced_lake & land
    # pass 2: lake cells add their open-water balance (often negative); flow is clipped at zero
    local2 = np.where(lake, gain_open, local)
    q = accumulate(order, rec, local2)
    marsh &= ~lake
    playa &= ~lake
    return {"runoff_mm": R, "discharge_m3s": q, "lake": lake, "terminal": terminal, "marsh": marsh, "playa": playa,
            "rec": rec,
            "filled": filled, "order": order, "lakes": lakes_info}


def stream_order(rec, order, q, threshold):
    """Strahler order of the channel network (cells with discharge >= threshold)."""
    n = rec.size
    so = np.zeros(n, np.int16)
    top = np.zeros(n, np.int16)       # highest order arriving
    cnt = np.zeros(n, np.int16)       # how many arrive with that order
    chan = q >= threshold
    for i in range(order.size - 1, -1, -1):
        c = order[i]
        if not chan[c]:
            continue
        if cnt[c] == 0:
            so[c] = 1
        else:
            so[c] = top[c] + (1 if cnt[c] >= 2 else 0)
        r = rec[c]
        if r >= 0 and chan[r]:
            if so[c] > top[r]:
                top[r] = so[c]
                cnt[r] = 1
            elif so[c] == top[r]:
                cnt[r] += 1
    return so

