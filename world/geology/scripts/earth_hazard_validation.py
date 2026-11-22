"""Test the seismic hazard index on Earth (experiment E16).

Sources: the plate boundary model PB2002 (Bird 2003; public copy github.com/fraxen/tectonicplates, folder
"original"): the boundary steps with their class, and the orogens (zones of distributed active deformation).
Reference: the hazard class of the GSHAP map (Giardini et al. 1999) at 71 cities (data/earth_hazard_sites.json).
Score: the share of cities whose modelled class is among the accepted ones.

Usage: python3 earth_hazard_validation.py --pb2002 DIR [--grid]   (DIR holds PB2002_steps.csv and
PB2002_orogens.dig.txt; --grid searches the weights and the decay length)
"""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath

import hazard

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"


RATE_BINS = (3.0, 10.0, 30.0)             # mm/yr: steps are grouped by class and relative velocity


def read_steps(d, step_km=10.0):
    """Points every ~step_km along the PB2002 boundary steps, by (class, velocity bin); and the length-weighted
    median velocity of each (class, bin)."""
    pts, vel = {}, {}
    lines = (Path(d) / "PB2002_steps.csv").read_text().splitlines()
    head = lines[0].split(",")
    ix = {k: head.index(k) for k in ("StartLong", "StartLat", "FinalLong", "FinalLat", "StepLength(km)", "StepClass",
                                     "VelocityLeft")}
    for ln in lines[1:]:
        f = ln.split(",")
        lo0, la0, lo1, la1 = (float(f[ix[k]]) for k in ("StartLong", "StartLat", "FinalLong", "FinalLat"))
        if lo1 - lo0 > 180:
            lo1 -= 360
        elif lo0 - lo1 > 180:
            lo1 += 360
        length = float(f[ix["StepLength(km)"]])
        v = float(f[ix["VelocityLeft"]])
        key = (f[ix["StepClass"]], int(np.digitize(v, RATE_BINS)))
        k = max(int(length / step_km), 1)
        t = np.linspace(0, 1, k + 1)
        p = pts.setdefault(key, ([], []))
        p[0].extend(la0 + t * (la1 - la0))
        p[1].extend(lo0 + t * (lo1 - lo0))
        vel.setdefault(key, []).append((v, length))
    med = {}
    for key, vl in vel.items():
        a = np.array(vl)
        o = np.argsort(a[:, 0])
        cw = np.cumsum(a[o, 1]) / a[:, 1].sum()
        med[key] = float(a[o, 0][np.searchsorted(cw, 0.5)])
    return pts, med


def read_orogens(d):
    """Orogen polygons [(lon, lat) array] of PB2002_orogens.dig.txt."""
    polys, cur = [], []
    for ln in (Path(d) / "PB2002_orogens.dig.txt").read_text().splitlines():
        s = ln.strip()
        if s.startswith("***"):
            if cur:
                polys.append(np.array(cur))
            cur = []
        elif s.startswith(("+", "-")) and "," in s:
            a, b = s.split(",")
            cur.append((float(a), float(b)))
    if cur:
        polys.append(np.array(cur))
    # longitudes made continuous: a polygon across the antimeridian would otherwise span the globe
    for p in polys:
        p[:, 0] = np.degrees(np.unwrap(np.radians(p[:, 0])))
    return polys


def inside(polys, lat, lon):
    m = np.zeros(len(lat), bool)
    for p in polys:
        q = p.copy()
        for shift in (0, 360, -360):
            m |= MPath(q + [shift, 0]).contains_points(np.c_[lon, lat])
    return m


def inputs(d, sites):
    lat = np.array([s[1] for s in sites])
    lon = np.array([s[2] for s in sites])
    pts, med = read_steps(d)
    sources = [(key[0], med[key], hazard.distance_km(np.array(v[0]), np.array(v[1]), lat, lon))
               for key, v in pts.items()]
    polys = read_orogens(d)
    zin = inside(polys, lat, lon)
    edge = np.concatenate([np.c_[np.linspace(p[i, 1], p[(i + 1) % len(p), 1], 6),
                                 np.linspace(p[i, 0], p[(i + 1) % len(p), 0], 6)]
                           for p in polys for i in range(len(p))])
    zdist = np.where(zin, 0.0, hazard.distance_km(edge[:, 0], edge[:, 1], lat, lon))
    return sources, zin, zdist


def score(cfg, sources, zin, zdist, sites):
    h = hazard.index(sources, zin * cfg["zone_weight"], zdist, cfg)
    c = hazard.classes(h, cfg)
    ok = np.array([ci in s[3] for ci, s in zip(c, sites)])
    return float(ok.mean()), c, h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb2002", required=True)
    ap.add_argument("--grid", action="store_true")
    args = ap.parse_args()
    cfg = json.loads((DATA / "geology.json").read_text())["hazard"]
    sites = json.loads((DATA / "earth_hazard_sites.json").read_text())["sites"]
    sources, zin, zdist = inputs(args.pb2002, sites)
    if args.grid:
        # boundary weights stay fixed by class (geology.json); searched: decay length, reference velocity, the
        # weight of the deforming zones and the three class thresholds
        grid = list(itertools.product((200, 300, 400, 500, 600), (10, 20), (0.45, 0.5), (0.2, 0.25, 0.3),
                                      (0.35, 0.4, 0.45), (0.55, 0.6, 0.65, 0.7)))
        ok = np.zeros((len(grid), len(sites)), bool)
        for g, (L, vr, wz, t1, t2, t3) in enumerate(grid):
            c2 = dict(cfg, decay_km=L, rate_ref_mm_yr=vr, zone_weight=wz, zone_weight_outside=wz,
                      class_thresholds=[t1, t2, t3])
            ok[g] = [ci in s[3] for ci, s in zip(score(c2, sources, zin, zdist, sites)[1], sites)]
        rate = ok.mean(1)
        for g in np.argsort(-rate, kind="stable")[:8]:
            print(round(float(rate[g]), 3), grid[g])
        # leave-one-out: choose the best setting without each city, score that city
        loo = [ok[np.argmax(np.delete(ok, i, 1).mean(1)), i] for i in range(len(sites))]
        print("leave-one-out hit rate", round(float(np.mean(loo)), 3))
        return
    hit, c, h = score(cfg, sources, zin, zdist, sites)
    print(json.dumps({"site_hit_rate": round(hit, 3), "sites": len(sites)}))
    for s, ci, hi in zip(sites, c, h):
        print(("  ok  " if ci in s[3] else "  MISS"), f"{s[0]:16s} index {hi:.2f} class {ci} accepted {s[3]}")


if __name__ == "__main__":
    main()
