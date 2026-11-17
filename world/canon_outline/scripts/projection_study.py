"""Experiment E02: how the flat map "The Known World" can lie on an Earth-sized sphere.

A flat map 11 800 x 7 800 km cannot keep scale, shape and north-up all at once on a sphere. This script measures the
trade-off for three ways of reading the map, and finds the latitudes at which canon climates fall into the latitude
bands where the same climates occur on Earth.

Ways of reading the map (inverse projections; the map's scale is fixed by the Wall in every case):
  A  cylindrical equidistant: north is up everywhere, north-south distances true, east-west scale true only at the
     Wall's latitude;
  M  Mercator, true at the Wall's latitude: north is up everywhere and shapes are kept (conformal), but the scale
     changes with latitude in every direction;
  B  equidistant conic, central meridian at the centre of the map: north-south distances true, east-west scale true
     at the Wall's latitude and at a second, optimised parallel; north tilts towards the map's edges;
  C  as B, but with the central meridian through the Wall (Westeros keeps north-up; Essos tilts more);
  S  sinusoidal, central meridian through the Wall: parallels horizontal, distances along every parallel and along
     the central meridian true, areas true; meridians lean away from the central meridian (shear);
  S2 as S, central meridian at the centre of the map.
For S the reported "scale" is the scale along the drawn vertical (h); east-west scale is exactly 1.

Input: the land mask from land_mask.py (white = land), 9050 x 6000 px.
Scale factor k = map length / true length (the map's scale being the Wall's). True area = map area / (h k).
Output: stats.json, scan.csv (one row per projection and Wall latitude), and a printed summary.
Deterministic: no randomness.
"""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

Image.MAX_IMAGE_PIXELS = None

R = 6371.0                       # km, Earth's mean radius (world/planet/PLANET.md)
KM_PER_PX = 1.3014               # E01: the Wall, 480 km over 368.8 px
WALL_WEST = (1568, 1015)
WALL_EAST = (1933, 962)
WALL_MID = ((WALL_WEST[0] + WALL_EAST[0]) / 2, (WALL_WEST[1] + WALL_EAST[1]) / 2)
WINTERFELL = (1550, 1607)        # map px, read on M1
STEP = 6                         # sample every STEP-th pixel

# Canon regions (boxes in map px, read on M1) and the latitude band where the same climate occurs on Earth.
# Sources of the bands: see experiments/E02_projection.md.
REGIONS = {
    #  name                 box (x0, y0, x1, y1)        band (deg N)   canon fact
    "haunted_forest":     ((1260, 730, 1910, 1000),   (50.0, 70.0)),   # boreal forest north of the Wall
    "dorne_desert":       ((1260, 3670, 2140, 3980),  (15.0, 35.0)),   # WE9, hot desert
    "red_waste":          ((5510, 3670, 6280, 4480),  (15.0, 35.0)),   # E6-E7, hot desert
    "dothraki_sea":       ((3830, 2680, 6120, 3830),  (30.0, 60.0)),   # E6, temperate grassland (steppe)
    "sothoryos_jungle":   ((4250, 5200, 5890, 5880),  (0.0, 15.0)),    # S2, tropical rainforest
    "summer_isles":       ((2300, 4820, 3180, 5740),  (0.0, 23.4)),    # I1, tropical
    "always_winter":      ((570, 115, 1760, 880),     (60.0, 90.0)),   # WE7, ice cap (Greenland 60-83 N)
}


def load_land(path):
    m = np.asarray(Image.open(path)) > 0
    if m.shape != (6000, 9050):
        raise SystemExit("expected the 9050 x 6000 land mask from land_mask.py")
    return m


def westeros_mask(land):
    lab, _ = ndimage.label(land)
    return lab == lab[WINTERFELL[1], WINTERFELL[0]]


def conic_params(phi1, phi2):
    if abs(phi2 - phi1) < 1e-9:
        n = math.sin(phi1)
    else:
        n = (math.cos(phi1) - math.cos(phi2)) / (phi2 - phi1)
    G = math.cos(phi1) / n + phi1
    return n, G


def inverse(kind, x, y, phi_wall, phi2=None, x0=None):
    """Map px -> (lat rad, lon rad, k = east-west scale factor, tilt of north in rad). h = 1 in all cases."""
    X = (x - (x0 if x0 is not None else WALL_MID[0])) * KM_PER_PX
    Y = (y - WALL_MID[1]) * KM_PER_PX            # downwards
    if kind == "A":
        lat = phi_wall - Y / R
        lon = X / (R * math.cos(phi_wall))
        k = math.cos(phi_wall) / np.cos(lat)
        return lat, lon, k, np.zeros_like(lat)
    if kind == "S":                              # sinusoidal: parallels horizontal and true, equal-area
        lat = phi_wall - Y / R
        lon = X / (R * np.cos(lat))
        tilt = np.arctan(lon * np.sin(lat))      # meridians lean away from the central meridian
        h = np.sqrt(1 + (lon * np.sin(lat)) ** 2)  # map length / true length along the drawn vertical
        return lat, lon, h, tilt
    if kind == "M":                              # Mercator, true scale at the Wall's latitude
        c = R * math.cos(phi_wall)
        yw = c * math.log(math.tan(math.pi / 4 + phi_wall / 2))
        lat = 2 * np.arctan(np.exp((yw - Y) / c)) - math.pi / 2
        lon = X / c
        k = math.cos(phi_wall) / np.cos(lat)     # same in every direction (conformal)
        return lat, lon, k, np.zeros_like(lat)
    n, G = conic_params(phi_wall, phi2)
    rho_w = R * (G - phi_wall)
    Xw = (WALL_MID[0] - x0) * KM_PER_PX
    if rho_w <= abs(Xw):
        raise ValueError("cone too steep: the Wall would lie beyond the apex circle")
    Ya = -math.sqrt(rho_w ** 2 - Xw ** 2)       # apex, above the Wall (Y measured from the Wall)
    rho = np.hypot(X, Y - Ya)
    lat = G - rho / R
    theta = np.arctan2(X, Y - Ya)
    lon = theta / n
    k = n * rho / (R * np.cos(lat))
    return lat, lon, k, theta


def best_phi2(kind, xs, ys, phi_wall, x0):
    best = None
    for p2 in np.radians(np.arange(5.0, 75.01, 0.5)):
        try:
            lat, _, k, _ = inverse(kind, xs, ys, phi_wall, p2, x0)
        except ValueError:
            continue
        if not np.all(np.isfinite(k)) or np.any(np.cos(lat) <= 0.02):
            continue
        err = float(np.sqrt(np.mean(np.log(k) ** 2)))
        if best is None or err < best[0]:
            best = (err, p2)
    return best[1]


def evaluate(kind, land, wes, phi_wall_deg, x0):
    ys, xs = np.nonzero(land[::STEP, ::STEP])
    ys = ys * STEP
    xs = xs * STEP
    phi_wall = math.radians(phi_wall_deg)
    phi2 = None if kind in ("A", "M", "S") else best_phi2(kind, xs, ys, phi_wall, x0)
    lat, lon, k, tilt = inverse(kind, xs, ys, phi_wall, phi2, x0)
    in_w = wes[ys, xs]
    res = {
        "projection": kind,
        "wall_lat": phi_wall_deg,
        "phi2": None if phi2 is None else round(math.degrees(phi2), 1),
        "scale_err_rms_all": round(float(np.sqrt(np.mean(np.log(k) ** 2))), 3),
        "scale_err_p95_all": round(float(np.percentile(np.abs(k - 1), 95)), 3),
        "scale_err_p95_westeros": round(float(np.percentile(np.abs(k[in_w] - 1), 95)), 3),
        "tilt_max_westeros_deg": round(float(np.degrees(np.abs(tilt[in_w]).max())), 1),
        "tilt_max_all_deg": round(float(np.degrees(np.abs(tilt).max())), 1),
        "lat_range_all": [round(float(np.degrees(lat.min())), 1), round(float(np.degrees(lat.max())), 1)],
        "lon_span_deg": round(float(np.degrees(lon.max() - lon.min())), 1),
        "westeros_area_mkm2": round(float(np.sum(1 / (k[in_w] * (k[in_w] if kind == "M" else 1)))
                                          * (STEP * KM_PER_PX) ** 2 / 1e6), 2) if kind != "S" else
                              round(float(in_w.sum() * (STEP * KM_PER_PX) ** 2 / 1e6), 2),
    }
    wf_lat, _, _, _ = inverse(kind, np.array([WINTERFELL[0]]), np.array([WINTERFELL[1]]), phi_wall, phi2, x0)
    res["winterfell_lat"] = round(float(np.degrees(wf_lat[0])), 1)
    fits = {}
    for name, ((bx0, by0, bx1, by1), (lo, hi)) in REGIONS.items():
        sel = (xs >= bx0) & (xs < bx1) & (ys >= by0) & (ys < by1)
        if not sel.any():
            fits[name] = None
            continue
        la = np.degrees(lat[sel])
        fits[name] = {
            "lat": [round(float(la.min()), 1), round(float(la.max()), 1)],
            "in_band": round(float(np.mean((la >= lo) & (la <= hi))), 2),
        }
    res["regions"] = fits
    res["min_in_band"] = min(v["in_band"] for v in fits.values() if v)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mask", required=True, help="land_mask.png from land_mask.py")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    land = load_land(args.mask)
    wes = westeros_mask(land)

    setups = {"A": None, "M": None, "B": 9050 / 2, "C": WALL_MID[0], "S": WALL_MID[0], "S2": 9050 / 2}
    rows = []
    best = {}
    for kind, x0 in setups.items():
        for wl in np.arange(50.0, 70.01, 0.5):
            r = evaluate(kind[0], land, wes, float(wl), x0)
            r["projection"] = kind
            rows.append(r)
            if kind not in best or r["min_in_band"] > best[kind]["min_in_band"]:
                best[kind] = r

    with open(out / "scan.csv", "w", newline="") as f:
        w = csv.writer(f)
        names = list(REGIONS)
        w.writerow(["projection", "wall_lat", "phi2", "scale_err_rms_all", "scale_err_p95_westeros",
                    "tilt_max_westeros_deg", "tilt_max_all_deg", "winterfell_lat", "min_in_band"] + names)
        for r in rows:
            w.writerow([r["projection"], r["wall_lat"], r["phi2"], r["scale_err_rms_all"],
                        r["scale_err_p95_westeros"], r["tilt_max_westeros_deg"], r["tilt_max_all_deg"],
                        r["winterfell_lat"], r["min_in_band"]] + [r["regions"][n]["in_band"] for n in names])
    (out / "stats.json").write_text(json.dumps(best, indent=2))
    print(json.dumps(best, indent=2))


if __name__ == "__main__":
    main()
