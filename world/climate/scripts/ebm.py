"""Seasonal two-dimensional energy balance model (EBM) on a latitude-longitude grid (experiment E09).

One temperature per cell (surface and lower atmosphere together, reduced to sea level):

    C dT/dt = S(lat, day) (1 - albedo) - (A + B T) + D lap(T) - C_atm u(lat) / (R cos lat) dT/dlon

- S: daily mean top-of-atmosphere insolation from the orbit (Berger 1978 formulae; Earth's orbit, world/planet).
- albedo: planetary albedo a0 + a2 P2(sin lat); snow/ice albedo where T < ice_T (land) or the sea is frozen.
- A + B T: outgoing longwave radiation (Budyko-type linear fit; North & Coakley 1979).
- D lap(T): meridional and zonal heat transport as diffusion on the unit sphere (North 1975), D in W m-2 K-1.
- advection of the atmospheric part (heat capacity C_atm) by a prescribed zonal wind u(lat): trade easterlies,
  mid-latitude westerlies, polar easterlies; this carries ocean temperatures onto the land downwind.
- C: ocean mixed layer (rho c h) or land (atmosphere + soil); frozen sea has the land value.

Backward Euler in time (dt = 1 day) with a sparse LU factorisation, rebuilt each month when the sea ice changes.
Deterministic.
"""

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

R_M = 6.371e6
DAY = 86400.0
YEAR_DAYS = 365
RHO_L = 917 * 3.34e5          # ice density x latent heat of fusion, J m-3
K_ICE = 2.2                   # thermal conductivity of sea ice, W m-1 K-1
MONTH_START = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])


def grid(res_deg):
    lat = 90 - (np.arange(int(180 / res_deg)) + 0.5) * res_deg
    lon = -180 + (np.arange(int(360 / res_deg)) + 0.5) * res_deg
    return lat, lon


def insolation(lat_deg, day, orbit):
    """Daily mean TOA insolation (W m-2) for each latitude on a given day of the year (day 0 = 1 January)."""
    e, eps, varpi = orbit["eccentricity"], math.radians(orbit["obliquity_deg"]), math.radians(orbit["perihelion_deg"])
    # mean longitude measured from the vernal equinox (day of equinox ~ 79), then true solar longitude
    L = 2 * math.pi * (day - orbit["equinox_day"]) / YEAR_DAYS
    lam = L + 2 * e * math.sin(L - varpi) + 1.25 * e * e * math.sin(2 * (L - varpi))
    dec = math.asin(math.sin(eps) * math.sin(lam))
    rho = ((1 + e * math.cos(lam - varpi)) / (1 - e * e)) ** 2
    phi = np.radians(lat_deg)
    x = np.clip(-np.tan(phi) * math.tan(dec), -1, 1)
    h0 = np.arccos(x)
    return orbit["solar_constant"] / math.pi * rho * (h0 * np.sin(phi) * math.sin(dec)
                                                     + np.cos(phi) * math.cos(dec) * np.sin(h0))


def zonal_wind(lat_deg, p):
    a = np.abs(lat_deg)
    return (p["westerly_ms"] * np.exp(-((a - p["westerly_lat"]) / p["westerly_width"]) ** 2)
            - p["trade_ms"] * np.exp(-((a - p["trade_lat"]) / p["trade_width"]) ** 2)
            - p["polar_ms"] * np.exp(-((a - p["polar_lat"]) / p["polar_width"]) ** 2))


def transport_operator(lat, lon, p):
    """Sparse matrix L with (L T) = D lap(T) - C_atm u / (R cos) dT/dlon, on a periodic lat-lon grid."""
    ny, nx = len(lat), len(lon)
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    phi = np.radians(lat)
    cos_c = np.cos(phi)
    cos_n = np.cos(phi + dphi / 2)       # face towards the north (lat decreasing index -> north is j-1)
    cos_s = np.cos(phi - dphi / 2)
    D = p["diffusion_W_m2_K"]
    u = zonal_wind(lat, p["wind"])
    idx = np.arange(ny * nx).reshape(ny, nx)
    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(r.ravel())
        cols.append(c.ravel())
        vals.append(np.broadcast_to(v, r.shape).ravel())

    # meridional diffusion: (1/cos) d/dphi (cos dT/dphi); row j-1 is to the north; no flux through the poles
    cn = (D * cos_n / (cos_c * dphi ** 2))[:, None] * np.ones((1, nx))
    cs_ = (D * cos_s / (cos_c * dphi ** 2))[:, None] * np.ones((1, nx))
    add(idx[1:], idx[:-1], cn[1:])
    add(idx[1:], idx[1:], -cn[1:])
    add(idx[:-1], idx[1:], cs_[:-1])
    add(idx[:-1], idx[:-1], -cs_[:-1])
    # zonal diffusion: (1/cos^2) d2T/dlam2, with cos limited near the poles (a polar filter)
    cz = np.maximum(cos_c, math.cos(math.radians(p["polar_cap_deg"])))
    kz = (p["zonal_diffusion_factor"] * D / (cz ** 2 * dlam ** 2))[:, None] * np.ones((1, nx))
    east = np.roll(idx, -1, 1)
    west = np.roll(idx, 1, 1)
    add(idx, east, kz)
    add(idx, west, kz)
    add(idx, idx, -2 * kz)
    # zonal advection, first-order upwind, of the atmospheric heat capacity
    ka = (p["C_atm"] * np.abs(u) / (R_M * cz * dlam))[:, None] * np.ones((1, nx))
    pos = (u >= 0)[:, None] * np.ones((1, nx), bool)
    add(idx[pos], west[pos], ka[pos])            # wind from the west: upstream is the western neighbour
    add(idx[~pos], east[~pos], ka[~pos])
    add(idx, idx, -ka)
    L = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(ny * nx,) * 2)
    return L


def gyre_streamfunction(lat, lon, land_frac, p):
    """Wind-driven ocean gyres (Stommel 1948): r lap(psi) + beta dpsi/dx = curl(tau) / rho0 on the ocean cells,
    psi = 0 on land. tau = rho_air Cd |u| u from the zonal wind profile. Returns psi (m3/s) on the grid."""
    g = p["ocean_gyre"]
    ny, nx = len(lat), len(lon)
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    phi = np.radians(lat)
    u = zonal_wind(lat, p["wind"])
    tau = g["rho_air"] * g["drag"] * np.abs(u) * u                        # N m-2, eastward
    # curl tau = -(1/(R cos)) d(tau cos)/dphi   (tau_y = 0); phi decreases with row index
    tc = tau * np.cos(phi)
    dtc = np.gradient(tc, -dphi)
    curl = -dtc / (R_M * np.maximum(np.cos(phi), 0.05))
    beta = 2 * 7.292e-5 * np.cos(phi) / R_M
    ocean = land_frac < 0.5
    idx = -np.ones((ny, nx), int)
    idx[ocean] = np.arange(ocean.sum())
    r = g["friction_per_s"]
    rows, cols, vals = [], [], []
    rhs = np.zeros(ocean.sum())
    cosc = np.maximum(np.cos(phi), math.cos(math.radians(p["polar_cap_deg"])))
    for j in range(ny):
        dx = R_M * cosc[j] * dlam
        dy = R_M * dphi
        for i in np.flatnonzero(ocean[j]):
            k = idx[j, i]
            diag = -2 * r / dx ** 2 - 2 * r / dy ** 2
            rows.append(k); cols.append(k); vals.append(diag)
            for jj, ii, coef in ((j, (i + 1) % nx, r / dx ** 2 + beta[j] / (2 * dx)),
                                 (j, (i - 1) % nx, r / dx ** 2 - beta[j] / (2 * dx)),
                                 (j - 1, i, r / dy ** 2), (j + 1, i, r / dy ** 2)):
                if 0 <= jj < ny and ocean[jj, ii]:
                    rows.append(k); cols.append(idx[jj, ii]); vals.append(coef)
            rhs[k] = curl[j] / g["rho_water"]
    A = sp.csr_matrix((vals, (rows, cols)), shape=(ocean.sum(),) * 2)
    psi = np.zeros((ny, nx))
    psi[ocean] = spla.spsolve(A.tocsc(), rhs)
    return psi


def gyre_operator(lat, lon, land_frac, p):
    """Sparse matrix for the heat carried by the gyres: -f rho c (U dT/dx + V dT/dy) on ocean cells (upwind),
    with U = -dpsi/dy, V = dpsi/dx (m2 s-1, depth-integrated) and f the share of the column that carries the
    surface temperature anomaly."""
    g = p["ocean_gyre"]
    ny, nx = len(lat), len(lon)
    psi = gyre_streamfunction(lat, lon, land_frac, p)
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    cosc = np.maximum(np.cos(np.radians(lat)), math.cos(math.radians(p["polar_cap_deg"])))
    dy = R_M * dphi
    dx = (R_M * cosc * dlam)[:, None]
    # y increases northward = decreasing row index
    U = -(np.roll(psi, 1, 0) - np.roll(psi, -1, 0)) / (2 * dy)
    U[0], U[-1] = 0, 0
    V = (np.roll(psi, -1, 1) - np.roll(psi, 1, 1)) / (2 * dx)
    w = (1 - land_frac) * g["heat_factor"] * 1025 * 3990
    idx = np.arange(ny * nx).reshape(ny, nx)
    east, west = np.roll(idx, -1, 1), np.roll(idx, 1, 1)
    north, south = np.roll(idx, 1, 0), np.roll(idx, -1, 0)
    ku = w * np.abs(U) / dx
    kv = w * np.abs(V) / dy
    kv[0], kv[-1] = 0, 0
    rows = [idx, idx, idx, idx]
    cols = [np.where(U >= 0, west, east), idx, np.where(V >= 0, south, north), idx]
    vals = [ku, -ku, kv, -kv]
    L = sp.csr_matrix((np.concatenate([v.ravel() for v in vals]),
                       (np.concatenate([r.ravel() for r in rows]), np.concatenate([c.ravel() for c in cols]))),
                      shape=(ny * nx,) * 2)
    return L, psi, U, V


def continentality(lat, lon, land_frac, p):
    """Moisture-source distance index M in [0, 1): 1 - exp(-d / L) on land, 0 at sea. d blends the distance to the
    sea measured upwind along the zonal wind (weight |u| / max|u|) with the straight distance to the sea."""
    ny, nx = len(lat), len(lon)
    res = abs(lat[1] - lat[0])
    dx = R_M / 1000 * np.radians(res) * np.maximum(np.cos(np.radians(lat)), 0.05)
    ocean = land_frac < 0.5
    u = zonal_wind(lat, p["wind"])
    d_up = np.full((ny, nx), np.inf)
    found = ocean.copy()
    d_up[ocean] = 0
    for k in range(1, nx):
        for j in range(ny):
            shift = k if u[j] >= 0 else -k          # wind from the west: look west (lower lon index)
            src = np.roll(ocean[j], shift)
            new = src & ~found[j]
            d_up[j][new] = k * dx[j]
            found[j] |= new
        if found.all():
            break
    # straight distance to the sea on the sphere
    from scipy.spatial import cKDTree
    LA, LO = np.meshgrid(np.radians(lat), np.radians(lon), indexing="ij")
    xyz = np.stack([np.cos(LA) * np.cos(LO), np.cos(LA) * np.sin(LO), np.sin(LA)], -1).reshape(-1, 3)
    d_iso, _ = cKDTree(xyz[ocean.ravel()]).query(xyz)
    d_iso = (d_iso * R_M / 1000).reshape(ny, nx)
    w = (np.abs(u) / np.abs(u).max())[:, None]
    d = w * np.minimum(d_up, 1e5) + (1 - w) * d_iso
    return np.where(ocean, 0.0, 1 - np.exp(-d / p["continentality"]["scale_km"]))


def qflux_field(lat, lon, land_frac, q):
    """Prescribed ocean heat-flux convergence ('Q-flux', as in slab-ocean models) for a warm poleward drift on the
    eastern side of one ocean basin: Q0 exp(-((lat - lat0)/width)^2) exp(-d_east / decay) on the basin's ocean cells,
    d_east = distance to the first land cell to the east along the parallel. The same amount of heat is taken
    uniformly from the basin's ocean between the equator and comp_lat_max, so the field sums to zero."""
    ny, nx = len(lat), len(lon)
    res = abs(lat[1] - lat[0])
    ocean = land_frac < 0.5
    lo = np.asarray(lon)
    if q["lon_min"] <= q["lon_max"]:
        inbox = (lo >= q["lon_min"]) & (lo <= q["lon_max"])
    else:                                            # the box wraps across 180 deg
        inbox = (lo >= q["lon_min"]) | (lo <= q["lon_max"])
    basin = ocean & inbox[None, :] & (lat[:, None] >= 0)
    dx = R_M / 1000 * np.radians(res) * np.maximum(np.cos(np.radians(lat)), 0.05)
    d_east = np.full((ny, nx), np.inf)
    land = ~ocean
    found = np.zeros((ny, nx), bool)
    for k in range(1, nx):
        hit = np.roll(land, -k, 1) & ~found          # land k cells to the east
        d_east[hit] = (k * dx[:, None] * np.ones((1, nx)))[hit]
        found |= hit
        if found.all():
            break
    Q = np.where(basin, q["Q0_W_m2"] * np.exp(-((lat[:, None] - q["lat0"]) / q["width_deg"]) ** 2)
                 * np.exp(-np.minimum(d_east, 1e5) / q["decay_km"]), 0.0)
    A = np.cos(np.radians(lat))[:, None] * np.ones((1, nx))
    comp = basin & (lat[:, None] <= q["comp_lat_max"])
    if comp.any():
        Q[comp] -= (Q * A).sum() / (A[comp]).sum()
    return Q


def run(land_frac, params, orbit, years=None, verbose=False, land_elev=None, ice_sheet=None):
    if params.get("two_layer"):
        return run_two_layer(land_frac, params, orbit, years, verbose, land_elev, ice_sheet)
    return run_one_layer(land_frac, params, orbit, years, verbose)


def run_one_layer(land_frac, params, orbit, years=None, verbose=False):
    """land_frac [ny, nx] on the grid of params['res_deg']. Returns monthly mean sea-level temperature (12, ny, nx)
    in deg C, the monthly sea-ice mask (12, ny, nx) and the daily insolation used (for checks)."""
    lat, lon = grid(params["res_deg"])
    ny, nx = len(lat), len(lon)
    n = ny * nx
    p = params
    L = transport_operator(lat, lon, p)
    if p.get("ocean_gyre", {}).get("heat_factor", 0) > 0:
        L = L + gyre_operator(lat, lon, land_frac, p)[0]
    P2 = 0.5 * (3 * np.sin(np.radians(lat)) ** 2 - 1)
    albedo_base = (p["albedo_a0"] + p["albedo_a2"] * P2)[:, None] * np.ones((1, nx))
    # far from moisture sources: fewer clouds (lower albedo) and drier air (weaker greenhouse, higher OLR)
    M = continentality(lat, lon, land_frac, p)
    albedo_base = albedo_base - p["continentality"]["albedo_drop"] * M
    A_eff = p["olr_A"] + p["continentality"]["olr_rise_W_m2"] * M
    C_ocean = 1025 * 3990 * p["mixed_layer_m"]
    C_land = p["C_land"]
    land = land_frac.astype(float)
    ice = np.zeros((ny, nx), bool)
    h_ice = np.zeros((ny, nx))
    years = int(years or p["spinup_years"])
    S_days = np.array([insolation(lat, d, orbit) for d in range(YEAR_DAYS)])        # (365, ny)
    # start from local radiative equilibrium of the annual mean (shortens the spin-up)
    T = np.clip((S_days.mean(0)[:, None] * (1 - albedo_base) - A_eff) / p["olr_B"], -40, 30) * np.ones((1, nx))
    monthly = np.zeros((12, ny, nx))
    ice_monthly = np.zeros((12, ny, nx), bool)
    dt = DAY
    for yr in range(years):
        for m in range(12):
            C = land * C_land + (1 - land) * np.where(ice, C_land * p["sea_ice_C_factor"], C_ocean)
            M = sp.diags((C / dt).ravel() + p["olr_B"]) - L
            lu = spla.splu(M.tocsc())
            acc = np.zeros((ny, nx))
            for d in range(MONTH_START[m], MONTH_START[m + 1]):
                snow = (land > 0.5) & (T < p["snow_T"])
                alb = np.where((h_ice > 0) | snow, np.maximum(albedo_base, p["ice_albedo"]), albedo_base)
                rhs = (C / dt) * T + S_days[d][:, None] * (1 - alb) - A_eff
                T = lu.solve(rhs.ravel()).reshape(ny, nx)
                # zero-layer sea ice (Semtner 1976): ice holds the surface at 0 C while it melts; open water
                # colder than freezing turns the deficit into ice; ice grows by conduction while the air is colder
                sea = land < 0.5
                melt = sea & (h_ice > 0) & (T > 0)
                h_ice[melt] -= C[melt] * T[melt] / RHO_L
                T[melt] = 0.0
                gone = melt & (h_ice < 0)
                T[gone] = -h_ice[gone] * RHO_L / C[gone]
                h_ice[gone] = 0.0
                grow = sea & (h_ice > 0) & (T < p["sea_ice_T"])
                h_ice[grow] += dt * K_ICE * (p["sea_ice_T"] - T[grow]) / (RHO_L * np.maximum(h_ice[grow], 0.1))
                form = sea & (h_ice <= 0) & (T < p["sea_ice_T"])
                h_ice[form] = C[form] * (p["sea_ice_T"] - T[form]) / RHO_L
                T[form] = p["sea_ice_T"]
                acc += T
            monthly[m] = acc / (MONTH_START[m + 1] - MONTH_START[m])
            ice = h_ice > 0
            ice_monthly[m] = ice
        if verbose:
            A = np.cos(np.radians(lat))[:, None] * np.ones((1, nx))
            print(yr, round(float((monthly.mean(0) * A).sum() / A.sum()), 2))
    return lat, lon, monthly, ice_monthly


def run_two_layer(land_frac, params, orbit, years=None, verbose=False, land_elev=None, ice_sheet=None):
    """Two layers per cell: surface Ts (land or mixed-layer ocean) and air Ta.

        Cs dTs/dt = (1 - f) S (1 - albedo) - g (Ts - Ta) - (a_w + b_w Ts)          surface
        Ca dTa/dt = f S (1 - albedo) + g (Ts - Ta) - (A + b_a Ta) + transport(Ta)    air

    f: share of absorbed sunlight taken up by the air; a_w + b_w Ts: surface emission escaping through the
    atmospheric window; g: surface-air exchange, large when the surface is warmer than the air (convection), small
    when it is colder (stable inversions, as over winter continents). Only the air is mixed and carried by the wind.

    Sea ice (zero-layer, Semtner 1976): where the sea is frozen, the surface temperature is that of the ice, from the
    balance of sunlight, exchange with the air, window emission and conduction through the ice from the water at
    freezing point below; it cannot exceed 0 C, and the surplus melts the ice. The water under the ice stays at
    freezing point.

    land_elev (ny, nx, m), optional: mean height of the land in each cell. Snow lies where the ground, not sea level,
    is below snow_T: Ts - lapse * land_elev (the snow line on high ground and ice-sheet plateaux).
    ice_sheet (ny, nx) in [0, 1], optional: share of the cell under an ice sheet, which keeps the ice albedo all year.
    Returns monthly mean surface temperature (deg C, sea level) and the sea-ice mask, as run_one_layer."""
    lat, lon = grid(params["res_deg"])
    ny, nx = len(lat), len(lon)
    n = ny * nx
    p = params
    tl = p["two_layer"]
    L = transport_operator(lat, lon, p)
    P2 = 0.5 * (3 * np.sin(np.radians(lat)) ** 2 - 1)
    albedo_base = (p["albedo_a0"] + p["albedo_a2"] * P2)[:, None] * np.ones((1, nx))
    M = continentality(lat, lon, land_frac, p)
    albedo_base = albedo_base - p["continentality"]["albedo_drop"] * M
    A_eff = p["olr_A"] + p["continentality"]["olr_rise_W_m2"] * M
    C_ocean = 1025 * 3990 * p["mixed_layer_m"]
    land = land_frac.astype(float)
    sea = land < 0.5
    Cs = land * p["C_land"] + (1 - land) * C_ocean
    Ca = p["C_atm"]
    f, aw, bw, ba = tl["air_absorption"], tl["window_a"], tl["window_b"], p["olr_B"]
    Tf = p["sea_ice_T"]
    qf = qflux_field(lat, lon, land_frac, p["qflux"]) if p.get("qflux", {}).get("Q0_W_m2", 0) else np.zeros((ny, nx))
    S_days = np.array([insolation(lat, d, orbit) for d in range(YEAR_DAYS)])
    Q0 = S_days.mean(0)[:, None] * (1 - albedo_base)
    Ta = np.clip((Q0 - A_eff - aw) / (ba + bw), -40, 30) * np.ones((1, nx))
    Ts = Ta.copy()
    z_snow = 0.0 if land_elev is None else p["snow_lapse_K_per_km"] / 1000 * land_elev
    sheet = np.zeros((ny, nx)) if ice_sheet is None else ice_sheet
    h_ice = np.zeros((ny, nx))
    years = int(years or p["spinup_years"])
    monthly = np.zeros((12, ny, nx))
    ice_monthly = np.zeros((12, ny, nx), bool)
    stable = np.zeros((ny, nx), bool)
    dt = DAY
    for yr in range(years):
        for m in range(12):
            g = np.where(stable, tl["exchange_stable"], tl["exchange_unstable"])
            G = sp.diags(g.ravel())
            top = sp.hstack([sp.diags((Cs / dt).ravel() + g.ravel() + bw), -G])
            bot = sp.hstack([-G, sp.diags(np.full(n, Ca / dt + ba)) + G - L])
            lu = spla.splu(sp.vstack([top, bot]).tocsc())
            acc = np.zeros((ny, nx))
            dsum = np.zeros((ny, nx))
            for d in range(MONTH_START[m], MONTH_START[m + 1]):
                ice = h_ice > 0
                snow = (land > 0.5) & (Ts - z_snow < p["snow_T"])
                alb = np.where(ice | snow, np.maximum(albedo_base, p["ice_albedo"]), albedo_base)
                alb = alb + sheet * (np.maximum(albedo_base, p["ice_albedo"]) - alb)
                Q = S_days[d][:, None] * (1 - alb)
                rhs = np.concatenate([((Cs / dt) * Ts + (1 - f) * Q - aw + qf).ravel(),
                                      ((Ca / dt) * Ta + f * Q - A_eff).ravel()])
                x = lu.solve(rhs)
                Ts = x[:n].reshape(ny, nx)
                Ta = x[n:].reshape(ny, nx)
                # sea ice: surface temperature from the ice's flux balance, growth by conduction, melt by surplus
                if ice.any():
                    kh = K_ICE / np.maximum(h_ice[ice], 0.05)
                    gi = g[ice]
                    fin = (1 - f) * Q[ice] + gi * Ta[ice] + kh * Tf - aw
                    Ti = fin / (gi + bw + kh)
                    surplus = np.where(Ti > 0, fin, 0.0)            # W m-2 left when the ice surface is at 0 C
                    Ti = np.minimum(Ti, 0.0)
                    hn = h_ice[ice] + dt * (kh * (Tf - Ti) - surplus - p["ocean_heat_under_ice_W_m2"]
                                            - np.maximum(qf[ice], 0)) / RHO_L
                    Ts[ice] = np.where(hn > 0, Ti, Tf)
                    h_ice[ice] = np.maximum(hn, 0.0)
                form = sea & (h_ice <= 0) & (Ts < Tf)
                h_ice[form] = C_ocean * (Tf - Ts[form]) / RHO_L
                Ts[form] = Tf
                acc += Ts
                dsum += Ts - Ta
            monthly[m] = acc / (MONTH_START[m + 1] - MONTH_START[m])
            stable = dsum < 0
            ice_monthly[m] = h_ice > 0
        if verbose:
            A = np.cos(np.radians(lat))[:, None] * np.ones((1, nx))
            print(yr, round(float((monthly.mean(0) * A).sum() / A.sum()), 2))
    return lat, lon, monthly, ice_monthly
