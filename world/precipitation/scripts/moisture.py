"""Monthly moisture-budget model for precipitation (experiment E10).

For each month the column water vapour W (kg m-2 = mm) is in steady state:

    div(W V) - div(K grad W) + k W = E,      P = k W

- V: low-level wind = zonal belts (trades, westerlies, polar easterlies) + meridional cells (Hadley, Ferrel, polar)
  about an ITCZ that moves with the season and towards the local thermal equator + a thermal component (towards
  warm ground and its geostrophic counterpart) + a Gill (1980) circulation forced by the tropical land-sea
  temperature contrast: the monsoons.
- K: eddy diffusivity (transient eddies carry moisture down its gradient; mid-latitude storm tracks).
- E: evaporation. Ocean: proportional to wind speed and saturation vapour density at the sea surface; none from
  sea ice. Land: potential evaporation limited by the water that fell (previous pass of the annual cycle).
- k = P / W: 1/tau0, larger in colder air, times exp(boost); boost from low-level convergence, upslope wind,
  minus the Hadley subsidence band, plus Sverdrup vertical motion (poleward flow rises, equatorward sinks).
  An optional dependence on column relative humidity (Bretherton et al. 2004) was tested and is off.

The finite-volume upwind divergence conserves water, so global precipitation equals global evaporation.
Deterministic.
"""

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

R_M = 6.371e6
DAY = 86400.0
MONTH_DAYS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])


def grid(res):
    lat = 90 - (np.arange(int(180 / res)) + 0.5) * res
    lon = -180 + (np.arange(int(360 / res)) + 0.5) * res
    return lat, lon


def qsat_rel(T):
    """Saturation vapour pressure relative to 20 C (Magnus formula)."""
    es = lambda t: 6.112 * np.exp(17.67 * t / (t + 243.5))
    return es(np.maximum(T, -60.0)) / es(20.0)


def box(a, n):
    """Repeated 3x3 box smoothing on a periodic-longitude grid."""
    for _ in range(n):
        p = np.pad(a, ((1, 1), (0, 0)), mode="edge")
        p = np.concatenate([p[:, -1:], p, p[:, :1]], 1)
        a = sum(p[1 + dj:p.shape[0] - 1 + dj, 1 + di:p.shape[1] - 1 + di] for dj in (-1, 0, 1) for di in (-1, 0, 1)) / 9
    return a


def winds(lat, lon, month, T, p, extra=None):
    """Low-level wind (u east, v north; m/s) for a month (0 = January) on the grid, T: surface temperature."""
    w = p["winds"]
    shift = w["itcz_mean_deg"] + w["itcz_amp_deg"] * math.cos(2 * math.pi * (month - 6.5) / 12)
    # the ITCZ follows the thermal equator, locally: far north over a hot summer continent (the South Asian and
    # West African monsoons), close to the equator over the oceans (Schneider, Bischoff & Haug 2014). Per
    # longitude: the temperature-weighted centroid latitude of the warmest air within 35 degrees of the equator,
    # smoothed along the circle; blended with the zonal-mean position by itcz_thermal_frac.
    trop = np.abs(lat) < 35
    Tt = box(T, w["thermal_smooth"])[trop]
    wt = np.exp((Tt - Tt.max(0)) / w["itcz_thermal_tau_K"])
    te = (wt * lat[trop][:, None]).sum(0) / wt.sum(0)
    k_lon = max(1, int(round(w["itcz_lon_smooth_deg"] / abs(lon[1] - lon[0]))))
    te = np.convolve(np.concatenate([te[-k_lon:], te, te[:k_lon]]), np.ones(2 * k_lon + 1) / (2 * k_lon + 1), "valid")
    f = w["itcz_thermal_frac"]
    x = lat[:, None] - ((1 - f) * shift + f * te)[None, :]
    a = np.abs(x)
    s = np.sign(x)
    # westerlies are stronger in the winter hemisphere (sharper equator-pole temperature contrast)
    season = math.cos(2 * math.pi * month / 12) * np.sign(lat)[:, None]     # +1 northern January, southern July
    wmult = 1 + w["westerly_winter_boost"] * season
    u = (w["westerly_ms"] * wmult * np.exp(-((a - 45) / 15) ** 2) - w["trade_ms"] * np.exp(-((a - 15) / 10) ** 2)
         - w["polar_ms"] * np.exp(-((a - 75) / 8) ** 2))
    v = np.where(a < 30, -w["hadley_ms"] * np.sin(np.pi * x / 30),
                 np.where(a < 60, w["ferrel_ms"] * s * np.sin(np.pi * (a - 30) / 30),
                          -w["polar_cell_ms"] * s * np.sin(np.pi * (np.minimum(a, 90) - 60) / 30)))
    # thermal wind towards warmer ground: k_th m/s per (K per 100 km), on smoothed temperature
    Ts = box(T, w["thermal_smooth"])
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    cosc = np.maximum(np.cos(np.radians(lat)), 0.1)[:, None]
    dTdy = (np.roll(Ts, 1, 0) - np.roll(Ts, -1, 0)) / (2 * R_M * dphi)
    dTdy[0], dTdy[-1] = 0, 0
    dTdx = (np.roll(Ts, -1, 1) - np.roll(Ts, 1, 1)) / (2 * R_M * cosc * dlam)
    u = u + w["thermal_k"] * dTdx * 1e5
    v = v + w["thermal_k"] * dTdy * 1e5
    # the same thermal pressure (low over warm ground, high over cool sea) in geostrophic balance: flow along the
    # isotherms, anticlockwise round warm lows in the north (clockwise in the south) - the summer monsoon and the
    # subtropical highs over cool oceans (dry west coasts, moist poleward flow on east coasts)
    s = (np.sign(lat) * np.minimum(np.abs(lat) / 10, 1))[:, None]
    u = u + w["geostrophic_k"] * s * dTdy * 1e5
    v = v - w["geostrophic_k"] * s * dTdx * 1e5
    # monsoon: the thermal pressure field of the land-sea contrast, p' = -alpha (T - zonal mean T) on smoothed
    # temperature (a heat low over a hot summer continent, a cold high over a winter one), drives a boundary-layer
    # wind in Ekman balance with linear drag r:  f k x V = -grad p'/rho - r V. Away from the equator it is nearly
    # geostrophic with cross-isobar inflow into the low; near the equator (f -> 0) it blows straight down the
    # pressure gradient, which carries the cross-equatorial monsoon flow (Somali jet) into the summer continent.
    if w.get("monsoon_alpha_pa_per_K", 0) > 0:
        Tm = box(T, w["monsoon_smooth"])
        pa = -w["monsoon_alpha_pa_per_K"] * (Tm - Tm.mean(1, keepdims=True))
        py = (np.roll(pa, 1, 0) - np.roll(pa, -1, 0)) / (2 * R_M * dphi)
        py[0], py[-1] = 0, 0
        px = (np.roll(pa, -1, 1) - np.roll(pa, 1, 1)) / (2 * R_M * cosc * dlam)
        f = (2 * 7.292e-5 * np.sin(np.radians(lat)))[:, None]
        r = w["monsoon_drag_per_s"]
        den = 1.2 * (r ** 2 + f ** 2)
        um = (-r * px - f * py) / den
        vm = (f * px - r * py) / den
        sp_m = np.hypot(um, vm)
        cap = np.minimum(1, w["monsoon_max_ms"] / np.maximum(sp_m, 1e-9))
        u = u + um * cap
        v = v + vm * cap
    if extra is not None:
        u = u + extra[0]
        v = v + extra[1]
    return u, v, x


def gill_operator(lat, lon, eps, c):
    """Sparse operator of the damped linear shallow-water (first baroclinic mode) problem on the sphere,
    Gill (1980) generalised from the equatorial beta-plane:
        eps u - f v + phi_x = 0,   eps v + f u + phi_y = 0,   eps phi + c^2 div(u, v) = -Q
    unknowns stacked [u, v, phi], centred differences on the grid, periodic in longitude, v = 0 beyond the poles."""
    ny, nx = len(lat), len(lon)
    n = ny * nx
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    phi = np.radians(lat)
    cos = np.maximum(np.cos(phi), 0.05)
    f = 2 * 7.292e-5 * np.sin(phi)
    idx = np.arange(n).reshape(ny, nx)
    E_ = np.roll(idx, -1, 1); W_ = np.roll(idx, 1, 1)
    N_ = np.vstack([idx[:1], idx[:-1]]); S_ = np.vstack([idx[1:], idx[-1:]])   # row 0 is the north
    gx = (1 / (2 * R_M * cos * dlam))[:, None] * np.ones((1, nx))
    gy = np.full((ny, nx), 1 / (2 * R_M * dphi))
    I = sp.identity(n, format="csr")
    def op(pos, neg, g, extra=None):
        r = np.concatenate([idx.ravel(), idx.ravel()])
        cidx = np.concatenate([pos.ravel(), neg.ravel()])
        v = np.concatenate([g.ravel(), -g.ravel()])
        return sp.csr_matrix((v, (r, cidx)), shape=(n, n))
    Dx = op(E_, W_, gx)
    Dy = op(N_, S_, gy)
    F = sp.diags(np.repeat(f, nx))
    # divergence with the metric term: (1/(R cos)) d(v cos)/dphi
    cosg = np.repeat(cos, nx)
    Dyc = sp.diags(1 / cosg) @ Dy @ sp.diags(cosg)
    A = sp.bmat([[eps * I, -F, Dx], [F, eps * I, Dy], [c ** 2 * Dx, c ** 2 * Dyc, eps * I]], format="csc")
    return spla.splu(A)


def gill_winds(lu, Q, ny, nx):
    n = ny * nx
    x = lu.solve(np.concatenate([np.zeros(2 * n), -Q.ravel()]))
    return x[:n].reshape(ny, nx), x[n:2 * n].reshape(ny, nx)


def divergence_operator(lat, lon, u, v):
    """Sparse D with (D W) = div(W V) (per second), first-order upwind finite volumes, no flux through the poles."""
    ny, nx = len(lat), len(lon)
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    phi = np.radians(lat)
    area = (R_M ** 2 * np.cos(phi) * dphi * dlam)[:, None] * np.ones((1, nx))
    idx = np.arange(ny * nx).reshape(ny, nx)
    rows, cols, vals = [], [], []

    def face(a, b, c, L):
        """flow c (m/s) across a face of length L from cell a to cell b (c > 0) or back (c < 0)."""
        fwd = c > 0
        src = np.where(fwd, a, b)
        dst = np.where(fwd, b, a)
        q = np.abs(c) * L
        rows.extend([src.ravel(), dst.ravel()])
        cols.extend([src.ravel(), src.ravel()])
        vals.extend([(q / area.ravel()[src.ravel()].reshape(q.shape)).ravel(),
                     (-q / area.ravel()[dst.ravel()].reshape(q.shape)).ravel()])

    # east faces: between i and i+1
    ue = 0.5 * (u + np.roll(u, -1, 1))
    face(idx, np.roll(idx, -1, 1), ue, np.full((ny, nx), R_M * dphi))
    # north faces: between row j (south) and row j-1 (north), for j = 1..ny-1
    vn = 0.5 * (v[1:] + v[:-1])
    Ln = (R_M * np.cos(phi[1:] + dphi / 2) * dlam)[:, None] * np.ones((1, nx))
    face(idx[1:], idx[:-1], vn, Ln)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(ny * nx,) * 2)


def diffusion_operator(lat, lon, K):
    """Sparse L with (L W) = div(K grad W) (per second), finite volumes, no flux through the poles: the transport of
    vapour by weather systems (eddies), down the gradient of W."""
    ny, nx = len(lat), len(lon)
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    phi = np.radians(lat)
    area = (R_M ** 2 * np.cos(phi) * dphi * dlam)[:, None] * np.ones((1, nx))
    idx = np.arange(ny * nx).reshape(ny, nx)
    rows, cols, vals = [], [], []

    def face(a, b, coef):
        for r_, c_, v_ in ((a, a, -coef / area.ravel()[a.ravel()].reshape(coef.shape)),
                           (a, b, coef / area.ravel()[a.ravel()].reshape(coef.shape)),
                           (b, b, -coef / area.ravel()[b.ravel()].reshape(coef.shape)),
                           (b, a, coef / area.ravel()[b.ravel()].reshape(coef.shape))):
            rows.append(r_.ravel()); cols.append(c_.ravel()); vals.append(v_.ravel())

    cosc = np.maximum(np.cos(phi), 0.1)[:, None] * np.ones((1, nx))
    # east faces: conductance K * L / dx with L = R dphi, dx = R cos dlam
    face(idx, np.roll(idx, -1, 1), K * (R_M * dphi) / (R_M * cosc * dlam))
    # north faces: L = R cos(phi_face) dlam, dy = R dphi
    Ln = (R_M * np.cos(phi[1:] + dphi / 2) * dlam)[:, None] * np.ones((1, nx))
    face(idx[1:], idx[:-1], K * Ln / (R_M * dphi))
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(ny * nx,) * 2)


def run(lat, lon, T_month, land_frac, elev, ice_month, p, passes=None):
    """T_month (12, ny, nx): sea-level surface temperature (C); elev: mean ground height (m, 0 at sea);
    ice_month (12, ny, nx) bool: sea ice. Returns P (12, ny, nx) mm/month, E likewise, and the winds."""
    ny, nx = len(lat), len(lon)
    pp = p["precip"]
    lapse = p["lapse_rate_K_per_km"] / 1000
    zs = box(elev, pp["orography_smooth"])
    dphi = math.radians(abs(lat[1] - lat[0]))
    dlam = math.radians(lon[1] - lon[0])
    cosc = np.maximum(np.cos(np.radians(lat)), 0.1)[:, None]
    dzdy = (np.roll(zs, 1, 0) - np.roll(zs, -1, 0)) / (2 * R_M * dphi)
    dzdy[0], dzdy[-1] = 0, 0
    dzdx = (np.roll(zs, -1, 1) - np.roll(zs, 1, 1)) / (2 * R_M * cosc * dlam)
    Ld = diffusion_operator(lat, lon, pp["eddy_diffusivity_m2_s"])
    passes = passes or pp["passes"]
    lu_g = gill_operator(lat, lon, 1 / (pp["gill_damping_days"] * DAY), pp["gill_c_ms"]) if pp.get("gill_q", 0) > 0 or pp.get("gill_t", 0) > 0 else None
    P = np.zeros((12, ny, nx))
    E = np.zeros((12, ny, nx))
    Pprev = np.full((12, ny, nx), 2.0)      # mm/day, first guess for land evaporation
    out_w = []
    diag = {"subsidence": np.zeros((ny, nx)), "sverdrup_w": np.zeros((ny, nx))}
    for ipass in range(passes):
        for m in range(12):
            Tg = T_month[m] - lapse * elev                  # ground-level temperature
            extra = None
            if lu_g is not None and pp.get("gill_t", 0) > 0:
                # heating of the land-sea temperature contrast (sensible heat of a summer continent, amplified by
                # the convection it organises) as the Gill forcing: Q ~ T - zonal mean T, on smoothed temperature
                Tm_ = box(T_month[m], 2)
                # Gill's theory is a tropical one: the forcing fades poleward of ~35 degrees
                Ta = (Tm_ - Tm_.mean(1, keepdims=True)) * np.exp(-(lat / pp["gill_lat_deg"]) ** 4)[:, None]
                ug, vg = gill_winds(lu_g, pp["gill_t"] * 1e-4 * Ta, ny, nx)
                spg = np.hypot(ug, vg)
                capg = np.minimum(1, pp["gill_max_ms"] / np.maximum(spg, 1e-9))
                extra = (ug * capg, vg * capg)
            elif lu_g is not None and ipass > 0:
                # latent heating of last pass's rain (anomaly from the zonal mean) drives a Gill circulation:
                # low-level inflow, westerlies on the west side of the heating (the south-west monsoon of the
                # Arabian Sea), easterlies on the east side
                Pa = Pprev[m] - Pprev[m].mean(1, keepdims=True)
                ug, vg = gill_winds(lu_g, pp["gill_q"] * 1e-4 * Pa, ny, nx)
                spg = np.hypot(ug, vg)
                capg = np.minimum(1, pp["gill_max_ms"] / np.maximum(spg, 1e-9))
                extra = (ug * capg, vg * capg)
            u, v, xi = winds(lat, lon, m, T_month[m], p, extra)
            speed = np.hypot(u, v)
            D = divergence_operator(lat, lon, u, v)
            conv = -(D @ np.ones(ny * nx)).reshape(ny, nx)   # -div(V) (1/s)
            w_oro = u * dzdx + v * dzdy                      # m/s, positive up the slope
            # the descending branch of the Hadley cell (|x| ~ 25 degrees from the local ITCZ): subsidence warms and
            # dries the free troposphere and caps convection under a trade inversion; the deserts of Earth lie
            # here. Stronger over land, where no sea-surface moisture feeds the boundary layer from below.
            sub = np.exp(-((np.abs(xi) - pp["subsidence_lat_deg"]) / pp["subsidence_width_deg"]) ** 2)
            sub = sub * (pp["subsidence_sea"] + (1 - pp["subsidence_sea"]) * land_frac)
            # Sverdrup vorticity balance in the lower troposphere, beta v = f dw/dz: poleward flow rises, equatorward
            # flow sinks (Rodwell & Hoskins 2001) - wet eastern and dry western flanks of the subtropical highs.
            # w at 1 km (m/s) = v H beta / f = v H / (R tan(lat)), latitude floored at 8 degrees
            tanl = np.tan(np.radians(np.maximum(np.abs(lat), 8)))[:, None]
            w_sv = v * np.sign(lat)[:, None] * 1000 / (R_M * tanl)
            boost = np.clip(pp["conv_k"] * conv * 1e6 + pp["oro_k"] * w_oro * 100 - pp["subsidence_k"] * sub
                            + pp["sverdrup_k"] * w_sv * 1000,
                            -pp["boost_max"], pp["boost_max"])
            k = (1 / (pp["tau0_days"] * DAY)) * (1 / qsat_rel(Tg)) ** pp["cold_beta"] * np.exp(boost)
            sea_e = np.where(ice_month[m], 0.0, pp["ocean_e0_mm_day"] * (speed + 1) / 6 * qsat_rel(T_month[m]))
            pet = pp["land_pet0_mm_day"] * qsat_rel(Tg) * (Tg > 0)
            land_e = np.minimum(pet, pp["land_recycling"] * Pprev[m])
            e = (land_frac * land_e + (1 - land_frac) * sea_e) / DAY        # mm/s
            # precipitation grows steeply with the column's relative humidity W / Wsat (Bretherton, Peters &
            # Back 2004: P ~ exp(15.6 (r - 0.603)) over tropical oceans); solved by fixed-point iteration
            Wsat = pp["wsat20_mm"] * qsat_rel(Tg)
            W = pp["rh_r0"] * Wsat
            M = D - Ld
            kk = k
            for _ in range(pp["rh_iterations"]):
                kk = k * np.exp(np.clip(pp["rh_a"] * (W / Wsat - pp["rh_r0"]), -pp["rh_clip"], pp["rh_clip"]))
                Wn = np.maximum(spla.spsolve((M + sp.diags(kk.ravel())).tocsc(), e.ravel()).reshape(ny, nx), 0)
                W = 0.5 * W + 0.5 * Wn
            # final, mass-conserving solve with the last precipitation rates (P = k W balances E exactly)
            k = kk
            W = np.maximum(spla.spsolve((M + sp.diags(k.ravel())).tocsc(), e.ravel()).reshape(ny, nx), 0)
            P[m] = k * W * DAY                                               # mm/day
            if ipass == passes - 1:
                diag["subsidence"] += sub / 12
                diag["sverdrup_w"] += w_sv / 12
            E[m] = e * DAY
            if len(out_w) < 12:
                out_w.append((u, v))
        Pprev = P.copy()
    return P * MONTH_DAYS[:, None, None], E * MONTH_DAYS[:, None, None], out_w, diag
