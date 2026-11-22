"""Surface mass balance and glacier geometry (experiment E11).

Positive-degree-day model, the standard of ice-sheet modelling (Reeh 1991; Huybrechts & de Wolde 1999; Calov &
Greve 2005):
- snowfall: the share of each month's precipitation that falls as snow, a linear ramp from all snow at
  snow_t0_C to all rain at snow_t1_C (monthly means; Marsiat 1994);
- melt: positive degree days of the month from its mean temperature and the standard deviation sigma of daily
  temperature about it (Calov & Greve 2005, closed form), melting snow first (ddf_snow) and then ice (ddf_ice);
- refreezing: up to refreeze_frac of the year's snowfall is retained as superimposed ice (Reeh 1991);
- elevation: a raised surface (ice) is colder by the lapse rate and drier by elevation desertification,
  precipitation x exp(gamma * dT) (Budd & Smith 1979; Huybrechts 2002, gamma ~ 0.07 per K).
All in metres of water equivalent per year. Deterministic.
"""

import math

import numpy as np
from scipy.special import erfc

SQ2PI = math.sqrt(2 * math.pi)


def pdd_month(T, sigma, days):
    """Expected positive degree days of a month with mean temperature T (C) (Calov & Greve 2005, eq. 6)."""
    s = sigma
    return days * (s / SQ2PI * np.exp(-T ** 2 / (2 * s ** 2)) + T / 2 * erfc(-T / (math.sqrt(2) * s)))


MONTH_DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])


def smb(T_month, P_month_mm, p, dz=None):
    """Annual surface mass balance (m w.e.) and its parts.
    T_month (12, ...) C at the reference surface, P_month_mm (12, ...) mm; dz (...) m: surface above the reference
    (ice thickness, or a sub-grid height), which lowers T by the lapse rate and P by elevation desertification."""
    shape = T_month.shape[1:]
    dT = 0 if dz is None else -p["lapse_K_per_km"] / 1000 * dz
    acc = np.zeros(shape)
    pdd = np.zeros(shape)
    for m in range(12):
        T = T_month[m] + dT
        P = P_month_mm[m] / 1000 * np.exp(p["desert_gamma_per_K"] * np.clip(dT, -30, 10))
        f = np.clip((p["snow_t1_C"] - T) / (p["snow_t1_C"] - p["snow_t0_C"]), 0, 1)
        acc += f * P
        pdd += pdd_month(T, p["sigma_K"], MONTH_DAYS[m])
    melt_snow_pot = p["ddf_snow_mm"] / 1000 * pdd
    melt_snow = np.minimum(melt_snow_pot, acc)
    pdd_left = pdd - melt_snow * 1000 / p["ddf_snow_mm"]
    melt_ice = p["ddf_ice_mm"] / 1000 * pdd_left
    refreeze = np.minimum(p["refreeze_frac"] * acc, melt_snow + melt_ice)
    b = acc - melt_snow - melt_ice + refreeze
    return b, {"accumulation": acc, "pdd": pdd}


def glacier_fraction(b0, b_up, dz_up, relief, aar):
    """Share of a cell covered by glaciers from its sub-grid relief.
    b0: balance at the cell mean height; b_up: balance dz_up higher. The balance is linear in height, the cell's
    heights spread uniformly over mean +- relief/2; glaciers cover the area above the equilibrium line divided by
    the accumulation-area ratio (glaciers reach below the ELA; AAR ~0.6, Meier & Post 1962)."""
    grad = (b_up - b0) / dz_up
    half = np.maximum(relief / 2, 1.0)
    lo = b0 - grad * half                                   # balance at the lowest and highest point of the cell
    hi = b0 + grad * half
    # share of the uniform height range where b > 0 (b is linear in height)
    with np.errstate(divide="ignore", invalid="ignore"):
        cross = np.clip(lo / (lo - hi), 0, 1)               # position of b = 0 from the bottom
    above = np.where((lo > 0) & (hi > 0), 1.0, np.where((lo <= 0) & (hi <= 0), 0.0,
                     np.where(hi > lo, 1 - cross, cross)))
    return np.clip(above / aar, 0, 1)
