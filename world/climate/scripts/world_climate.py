"""The world's climate models in one place: temperature (E09) and precipitation (E10) on the cubed sphere.

Used by the temperature layer (climate.py), the precipitation layer (precipitation.py) and the ice layer
(world/ice/scripts/ice.py), so that all three see exactly the same climate for the same surface.

The surface is given as the ground height of every 10 km cell (bedrock plus ice, 0 at sea) and, optionally, the
share of each cell under an ice sheet (bright all year). Deterministic.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
sys.path.insert(0, str(WORLD / "precipitation" / "scripts"))
import cubesphere as cs     # noqa: E402
import ebm                  # noqa: E402
import moisture             # noqa: E402

TCFG = WORLD / "climate" / "data" / "climate.json"
PCFG = WORLD / "precipitation" / "data" / "precipitation.json"


def cell_latlon_deg(n=cs.N):
    clat = np.zeros((6, n, n))
    clon = np.zeros((6, n, n))
    for f in range(6):
        la, lo = cs.cell_latlon(f, n)
        clat[f], clon[f] = np.degrees(la), np.degrees(lo)
    return clat, clon


def coarse_mean(field_cube, lat, lon, sub=4):
    """Mean of a cube field over each cell of a coarse latitude-longitude grid (sub x sub samples per cell)."""
    res = lat[0] - lat[1]
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros((len(lat), len(lon)))
    for a in off:
        for b in off:
            LA, LO = np.meshgrid(np.clip(lat + a * res, -89.99, 89.99), lon + b * res, indexing="ij")
            f, j, i = cs.latlon_to_cell(np.radians(LA), np.radians(LO), cs.N)
            acc += field_cube[f, j, i]
    return acc / sub ** 2


def to_cube(field, lat, lon, cell_lat, cell_lon):
    """Bilinear interpolation of a coarse [ny, nx] field (lat descending, periodic lon) to cell centres."""
    res = lat[0] - lat[1]
    y = np.clip((lat[0] - cell_lat) / res, 0, len(lat) - 1.000001)
    x = (cell_lon - lon[0]) / res
    j0 = np.floor(y).astype(int)
    fy = y - j0
    j1 = np.minimum(j0 + 1, len(lat) - 1)
    i0 = np.floor(x).astype(int) % len(lon)
    fx = x - np.floor(x)
    i1 = (i0 + 1) % len(lon)
    return ((1 - fy) * ((1 - fx) * field[j0, i0] + fx * field[j0, i1])
            + fy * ((1 - fx) * field[j1, i0] + fx * field[j1, i1]))


def temperature(ground, land, ice_sheet=None, clat=None, clon=None):
    """Monthly temperature. ground: height of the land surface (m, 0 at sea), land: bool cube, ice_sheet: share of
    each cube cell under an ice sheet or None. Returns a dict: lat, lon, sea-level monthly field and sea ice on the
    coarse grid, and t_cube (12, 6, n, n) float32, ground-level temperature in C."""
    cfg = json.loads(TCFG.read_text())
    p = cfg["ebm"]
    # the warm drift (Q-flux) is calibrated on the North Atlantic; on the world it lies in the Sunset Sea
    p["qflux"] = dict(p["qflux"], lon_min=cfg["world_qflux_box"]["lon_min"], lon_max=cfg["world_qflux_box"]["lon_max"])
    if clat is None:
        clat, clon = cell_latlon_deg()
    lat, lon = ebm.grid(p["res_deg"])
    lf = coarse_mean(land.astype(float), lat, lon)
    zsum = coarse_mean(np.where(land, ground, 0.0), lat, lon)
    zl = np.where(lf > 0.01, zsum / np.maximum(lf, 0.01), 0.0)          # mean height of the land in the cell
    sheet = None if ice_sheet is None else coarse_mean(np.where(land, ice_sheet, 0.0), lat, lon)
    lat, lon, monthly, seaice = ebm.run(lf, p, cfg["orbit"], land_elev=zl, ice_sheet=sheet)
    lapse = cfg["lapse_rate_K_per_km"] / 1000
    g = np.where(land, ground, 0.0)
    # statistical correction fitted on Earth's stations (earth_validation.py, fit_correction): on land subtract
    # b0 + b1 z_km + b2 M, with the coefficients moving from winter to summer through the year
    sc = cfg["station_correction"]
    M = to_cube(ebm.continentality(lat, lon, lf, p), lat, lon, clat, clon)
    zk = np.minimum(g / 1000, sc["z_cap_km"])
    bw, bs = np.array(sc["winter"]), np.array(sc["summer"])
    hemi = np.sign(clat) * np.minimum(np.abs(clat) / sc["equator_ramp_deg"], 1)
    t = np.zeros((12,) + land.shape, np.float32)
    for m in range(12):
        s = np.cos(2 * np.pi * m / 12) * hemi                  # +1: northern January / southern July (winter)
        b = [(bw[k] + bs[k]) / 2 + (bw[k] - bs[k]) / 2 * s for k in range(3)]
        corr = np.where(land, b[0] + b[1] * zk + b[2] * M, 0.0)
        t[m] = to_cube(monthly[m], lat, lon, clat, clon) - lapse * g - corr
    return {"lat": lat, "lon": lon, "monthly": monthly, "sea_ice": seaice, "land_frac": lf, "t_cube": t}


def precipitation(temp, ground, land, clat=None, clon=None):
    """Monthly precipitation from the moisture model, driven by temperature(...)'s coarse fields.
    Returns a dict: lat, lon, P (12, ny, nx) mm on the model grid, and p_cube (12, 6, n, n) float32 mm."""
    cfg = json.loads(PCFG.read_text())
    if clat is None:
        clat, clon = cell_latlon_deg()
    lat2, lon2, T2, ice2 = temp["lat"], temp["lon"], temp["monthly"], temp["sea_ice"]
    res = cfg["precip"]["res_deg"]
    lat, lon = moisture.grid(res)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    T = np.stack([to_cube(T2[m], lat2, lon2, LA, LO) for m in range(12)])
    ice = np.stack([to_cube(ice2[m].astype(float), lat2, lon2, LA, LO) > 0.5 for m in range(12)])
    lf = coarse_mean(land.astype(float), lat, lon)
    z = coarse_mean(np.where(land, ground, 0.0), lat, lon)
    P, E, _, diag = moisture.run(lat, lon, T, lf, z, ice, cfg)
    pc = np.zeros((12,) + land.shape, np.float32)
    for m in range(12):
        pc[m] = np.maximum(to_cube(P[m], lat, lon, clat, clon), 0)
    return {"lat": lat, "lon": lon, "P": P, "p_cube": pc}
