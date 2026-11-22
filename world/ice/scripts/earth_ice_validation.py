"""Test the mass-balance model on Earth: where does present-day climate keep ice? (experiment E11)

Earth inputs, all made the same way as for the world:
- temperature: the E09 energy balance model on Earth's land mask, sea level, 2 degrees; to 10' bilinearly, minus
  the lapse rate times ETOPO1 surface height, minus the E09 station correction;
- precipitation: the E10 moisture model on Earth, 2 degrees, to 10' bilinearly;
- surface: ETOPO1 ice surface at 10' (the present ice sheets' own surface, so this tests the balance, not the
  ice thickness); sub-grid relief from the spread of heights in the 3 x 3 neighbourhood;
- observed ice: Natural Earth 1:10m glaciated areas (public domain), rasterised at 10'.
Scores per region: observed and modelled ice area, and over all land the hit rate and false-alarm area.

Usage: python3 earth_ice_validation.py --earth-data DIR [--out DIR] [--set smb.key=value ...]
"""

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MPath
from scipy.interpolate import RegularGridInterpolator

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
DATA = HERE.parent / "data"
sys.path.insert(0, str(WORLD / "climate" / "scripts"))
sys.path.insert(0, str(WORLD / "precipitation" / "scripts"))
import ebm                                   # noqa: E402
import earth_validation as tval              # noqa: E402
import earth_precip_validation as pval       # noqa: E402
import smb as smbm                           # noqa: E402

NE_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/"
          "ne_10m_glaciated_areas.geojson")
RES = 1 / 6
REGIONS = {   # name: (lat_min, lat_max, lon_min, lon_max)
    "Antarctica": (-90, -60, -180, 180),
    "Greenland": (59, 84, -75, -10),
    "Canadian Arctic islands": (72, 84, -125, -60),
    "Baffin Island": (62, 74, -82, -60),
    "Iceland": (63, 67, -25, -13),
    "Svalbard and Franz Josef": (76, 82, 10, 65),
    "Novaya Zemlya and Severnaya Zemlya": (70, 82, 50, 110),
    "Alaska and Yukon": (55, 66, -160, -128),
    "Patagonia": (-56, -40, -76, -68),
    "Scandinavia": (58, 71, 4, 30),
    "Alps": (44, 48.5, 5, 16),
    "High Mountain Asia": (27, 46, 65, 105),
    "Siberia (no ice)": (55, 72, 60, 140),
    "Canadian mainland (no ice)": (50, 68, -120, -65),
}


def grid10():
    lat = 90 - (np.arange(1080) + 0.5) * RES
    lon = -180 + (np.arange(2160) + 0.5) * RES
    return lat, lon


def etopo10(earth_data):
    import h5py
    f = Path(earth_data) / "earth-topography-10arcmin.nc"
    if not f.exists():
        Path(earth_data).mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(pval.TOPO_URL, f)
    z = h5py.File(f, "r")["topography"][:].astype(float)            # nodes, south first, 1081 x 2161
    zc = 0.25 * (z[:-1, :-1] + z[1:, :-1] + z[:-1, 1:] + z[1:, 1:])   # cell centres
    return zc[::-1]


def observed_ice(earth_data, lat, lon):
    f = Path(earth_data) / "ne_10m_glaciated_areas.geojson"
    if not f.exists():
        urllib.request.urlretrieve(NE_URL, f)
    gj = json.loads(f.read_text())
    mask = np.zeros((len(lat), len(lon)), bool)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    for feat in gj["features"]:
        g = feat["geometry"]
        polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:
            ring = np.array(poly[0])
            lo0, la0 = ring.min(0)
            lo1, la1 = ring.max(0)
            js = np.where((lat >= la0 - RES) & (lat <= la1 + RES))[0]
            is_ = np.where((lon >= lo0 - RES) & (lon <= lo1 + RES))[0]
            if len(js) == 0 or len(is_) == 0:
                continue
            sub = np.c_[LO[np.ix_(js, is_)].ravel(), LA[np.ix_(js, is_)].ravel()]
            inside = MPath(ring).contains_points(sub).reshape(len(js), len(is_))
            for hole in poly[1:]:
                inside &= ~MPath(np.array(hole)).contains_points(sub).reshape(len(js), len(is_))
            mask[np.ix_(js, is_)] |= inside
    return mask


def to10(field, lat_c, lon_c, lat, lon):
    """Bilinear from the coarse grid (lat descending, periodic lon) to the 10' grid."""
    fp = np.concatenate([field[:, -1:], field, field[:, :1]], 1)
    lp = np.concatenate([[lon_c[0] - (lon_c[1] - lon_c[0])], lon_c, [lon_c[-1] + (lon_c[1] - lon_c[0])]])
    it = RegularGridInterpolator((lat_c[::-1], lp), fp[::-1], bounds_error=False, fill_value=None)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    return it((np.clip(LA, lat_c[-1], lat_c[0]), LO))


def load_earth(earth_data):
    """Earth climate on the 10' grid (expensive: the temperature and moisture models; cache the result)."""
    tcfg = json.loads((WORLD / "climate" / "data" / "climate.json").read_text())
    pcfg = json.loads((WORLD / "precipitation" / "data" / "precipitation.json").read_text())
    lat, lon = grid10()
    z = etopo10(earth_data)
    land = tval.globe.is_land(*np.meshgrid(lat, lon, indexing="ij")) | (lat[:, None] < -60) & (z > 0)
    obs = observed_ice(earth_data, lat, lon) & land
    res2 = tcfg["ebm"]["res_deg"]
    k = int(round(res2 / RES))
    lf2 = tval.earth_land_fraction(res2)
    zland = np.where(land, np.maximum(z, 0), 0).reshape(len(lat) // k, k, len(lon) // k, k).sum((1, 3))
    nland = land.reshape(len(lat) // k, k, len(lon) // k, k).sum((1, 3))
    zl2 = np.where(nland > 0, zland / np.maximum(nland, 1), 0)
    sheet2 = obs.reshape(len(lat) // k, k, len(lon) // k, k).mean((1, 3))
    # the temperature model sees the land height (snow line) and the observed ice sheets (albedo), as on the world
    lat2, lon2, Tm, ice = ebm.run(lf2, tcfg["ebm"], tcfg["orbit"], land_elev=zl2, ice_sheet=sheet2)
    M2 = ebm.continentality(lat2, lon2, lf2, tcfg["ebm"])
    latp, lonp = pval.moisture.grid(pcfg["precip"]["res_deg"])
    pe = {"T": np.stack([pval.to_grid(Tm[m], lat2, lon2, latp, lonp) for m in range(12)]),
          "ice": np.stack([pval.to_grid(ice[m].astype(float), lat2, lon2, latp, lonp) > 0.5 for m in range(12)]),
          "land": tval.earth_land_fraction(pcfg["precip"]["res_deg"]),
          "elev": pval.earth_topography(Path(earth_data), pcfg["precip"]["res_deg"])}
    res, P, E = pval.evaluate(pcfg, pe)
    ground = np.where(land, np.maximum(z, 0), 0)
    sc = tcfg["station_correction"]
    M = to10(M2, lat2, lon2, lat, lon)
    zk = np.minimum(ground / 1000, sc["z_cap_km"])
    bw, bs = np.array(sc["winter"]), np.array(sc["summer"])
    hemi = (np.sign(lat) * np.minimum(np.abs(lat) / sc["equator_ramp_deg"], 1))[:, None]
    lapse = tcfg["lapse_rate_K_per_km"] / 1000
    T = np.zeros((12, len(lat), len(lon)), np.float32)
    Pm = np.zeros((12, len(lat), len(lon)), np.float32)
    for m in range(12):
        s = np.cos(2 * np.pi * m / 12) * hemi
        b = [(bw[k] + bs[k]) / 2 + (bw[k] - bs[k]) / 2 * s for k in range(3)]
        corr = np.where(land, b[0] + b[1] * zk + b[2] * M, 0.0)
        T[m] = to10(Tm[m], lat2, lon2, lat, lon) - lapse * ground - corr
        Pm[m] = np.maximum(to10(P[m], latp, lonp, lat, lon), 0)
    # sub-grid relief: spread of the 3 x 3 neighbourhood of 10' heights
    zz = np.pad(ground, ((1, 1), (0, 0)), mode="edge")
    stack = [np.roll(zz, d, 1)[1 + a:1 + a + len(lat)] for a in (-1, 0, 1) for d in (-1, 0, 1)]
    relief = np.max(stack, 0) - np.min(stack, 0)
    return {"lat": lat, "lon": lon, "z": ground, "land": land, "T": T, "P": Pm, "relief": relief,
            "obs": obs}


def evaluate(cfg, earth):
    lat, lon = earth["lat"], earth["lon"]
    land = earth["land"]
    p, g = cfg["smb"], cfg["glacier"]
    b0, _ = smbm.smb(earth["T"], earth["P"], p)
    bu, _ = smbm.smb(earth["T"], earth["P"], p, dz=np.full(b0.shape, float(g["dz_up_m"])))
    frac = np.where(land, smbm.glacier_fraction(b0, bu, g["dz_up_m"], earth["relief"], g["aar"]), 0)
    A = (np.cos(np.radians(lat)) * (111.2 * RES) ** 2)[:, None] * np.ones((1, len(lon)))
    obs = earth["obs"]
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    rows = {}
    for name, (a0, a1, o0, o1) in REGIONS.items():
        m = land & (LA >= a0) & (LA <= a1) & (LO >= o0) & (LO <= o1)
        if name == "Antarctica":
            m = land & (LA <= a1)
        rows[name] = [round(float((A * obs)[m].sum()) / 1e3), round(float((A * frac)[m].sum()) / 1e3)]
    model = frac > 0.5
    res = {
        "observed_ice_Mkm2": round(float((A * obs)[land].sum()) / 1e6, 2),
        "model_ice_Mkm2": round(float((A * frac)[land].sum()) / 1e6, 2),
        "hit_rate": round(float((A * (obs & model)).sum() / (A * obs).sum()), 3),
        "false_alarm_Mkm2": round(float((A * (model & ~obs & land)).sum()) / 1e6, 2),
        "regions_kkm2_obs_model": rows,
    }
    return res, b0, frac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--earth-data", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--set", nargs="*", default=[])
    args = ap.parse_args()
    cfg = json.loads((DATA / "ice.json").read_text())
    for kv in args.set:
        k, v = kv.split("=")
        a, q = k.split(".")
        cfg[a][q] = float(v)
    earth = load_earth(args.earth_data)
    res, b0, frac = evaluate(cfg, earth)
    print(json.dumps(res, indent=1))
    if args.out:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "earth_ice_validation.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()


# ---------------- the coupled ice-climate equilibrium on Earth (hysteresis test) ----------------

def latlon_neighbours(ny, nx):
    """4-neighbour flat indices on the 10' grid (periodic in longitude; at the poles a cell is its own neighbour)."""
    idx = np.arange(ny * nx).reshape(ny, nx)
    e = np.roll(idx, -1, 1)
    w = np.roll(idx, 1, 1)
    nn = np.vstack([idx[:1], idx[:-1]])
    ss = np.vstack([idx[1:], idx[-1:]])
    return np.stack([e, w, ss, nn]).reshape(4, -1).astype(np.int32)


def earth_geometry(earth_data):
    """Earth on the 10' grid for ice.equilibrate: unloaded bed = ETOPO1 surface off the ice and a uniform
    bed_under_ice_m under today's ice (Earth's bed under the sheets is not in the inputs; the rebounded beds of
    Greenland and East Antarctica are a few hundred metres above sea level), and a climate function that reruns the
    temperature and moisture models for any surface and ice cover."""
    cfg = json.loads((DATA / "ice.json").read_text())
    tcfg = json.loads((WORLD / "climate" / "data" / "climate.json").read_text())
    pcfg = json.loads((WORLD / "precipitation" / "data" / "precipitation.json").read_text())
    lat, lon = grid10()
    ny, nx = len(lat), len(lon)
    z = etopo10(earth_data)
    land = tval.globe.is_land(*np.meshgrid(lat, lon, indexing="ij")) | (lat[:, None] < -60) & (z > 0)
    obs = observed_ice(earth_data, lat, lon) & land
    bed = np.where(land, np.where(obs, cfg["earth_test"]["bed_under_ice_m"], np.maximum(z, 0)), 0.0)
    res2 = tcfg["ebm"]["res_deg"]
    k = int(round(res2 / RES))
    lf2 = tval.earth_land_fraction(res2)
    latp, lonp = pval.moisture.grid(pcfg["precip"]["res_deg"])
    lfp = tval.earth_land_fraction(pcfg["precip"]["res_deg"])
    kp = int(round(pcfg["precip"]["res_deg"] / RES))
    nland = land.reshape(ny // k, k, nx // k, k).sum((1, 3))
    sc = tcfg["station_correction"]
    bw, bs = np.array(sc["winter"]), np.array(sc["summer"])
    hemi = (np.sign(lat) * np.minimum(np.abs(lat) / sc["equator_ramp_deg"], 1))[:, None]
    lapse = tcfg["lapse_rate_K_per_km"] / 1000
    zz = np.pad(np.where(land, np.maximum(z, 0), 0), ((1, 1), (0, 0)), mode="edge")
    stack = [np.roll(zz, d, 1)[1 + a:1 + a + ny] for a in (-1, 0, 1) for d in (-1, 0, 1)]
    relief = np.where(obs, 0.0, np.max(stack, 0) - np.min(stack, 0))

    def climate(surface, share):
        g = np.where(land, surface, 0.0)
        zl2 = np.where(nland > 0, g.reshape(ny // k, k, nx // k, k).sum((1, 3)) / np.maximum(nland, 1), 0)
        sh2 = None if share is None else np.where(land, share, 0).reshape(ny // k, k, nx // k, k).mean((1, 3))
        lat2, lon2, Tm, ice = ebm.run(lf2, tcfg["ebm"], tcfg["orbit"], land_elev=zl2, ice_sheet=sh2)
        M = to10(ebm.continentality(lat2, lon2, lf2, tcfg["ebm"]), lat2, lon2, lat, lon)
        pe = {"T": np.stack([pval.to_grid(Tm[m], lat2, lon2, latp, lonp) for m in range(12)]),
              "ice": np.stack([pval.to_grid(ice[m].astype(float), lat2, lon2, latp, lonp) > 0.5 for m in range(12)]),
              "land": lfp, "elev": g.reshape(ny // kp, kp, nx // kp, kp).mean((1, 3))}
        P, E, _, _ = pval.moisture.run(latp, lonp, pe["T"], pe["land"], pe["elev"], pe["ice"], pcfg)
        zk = np.minimum(g / 1000, sc["z_cap_km"])
        T = np.zeros((12, ny, nx), np.float32)
        Pm = np.zeros((12, ny, nx), np.float32)
        for m in range(12):
            s = np.cos(2 * np.pi * m / 12) * hemi
            b = [(bw[q] + bs[q]) / 2 + (bw[q] - bs[q]) / 2 * s for q in range(3)]
            T[m] = to10(Tm[m], lat2, lon2, lat, lon) - lapse * g - np.where(land, b[0] + b[1] * zk + b[2] * M, 0)
            Pm[m] = np.maximum(to10(P[m], latp, lonp, lat, lon), 0)
        return T, Pm

    A = (np.cos(np.radians(lat)) * (111.2 * RES) ** 2)[:, None] * np.ones((1, nx))
    dlon = (np.cos(np.radians(lat)) * 111.2e3 * RES)[:, None] * np.ones((1, nx))
    dlat = np.full((ny, nx), 111.2e3 * RES)
    step = np.stack([dlon, dlon, dlat, dlat]).reshape(4, -1)
    return {"land": land, "bed": bed, "relief": relief, "A": A, "step": step,
            "sea_depth": np.where(land, 0.0, np.maximum(-z, 0)),
            "nb": latlon_neighbours(ny, nx), "climate": climate, "obs": obs, "lat": lat, "lon": lon}


def hysteresis(earth_data, start, cfg=None):
    import ice as icem
    cfg = cfg or json.loads((DATA / "ice.json").read_text())
    geom = earth_geometry(earth_data)
    res = icem.equilibrate(geom, cfg, start)
    lat, lon, A, land, obs = geom["lat"], geom["lon"], geom["A"], geom["land"], geom["obs"]
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    rows = {}
    for name, (a0, a1, o0, o1) in REGIONS.items():
        m = land & (LA >= a0) & (LA <= a1) & (LO >= o0) & (LO <= o1)
        rows[name] = [round(float(A[m & obs].sum()) / 1e3), round(float(A[m & res["sheet"]].sum()) / 1e3)]
    out = {"start": start, "history": res["history"], "regions_kkm2_obs_sheet": rows,
           "max_thickness_m": int(res["H"].max()),
           "greenland_max_surface_m": int(res["surface"][(LA > 59) & (LO > -75) & (LO < -10)].max()),
           "antarctica_mean_surface_m": int(res["surface"][(LA < -60) & res["sheet"]].mean())
           if (res["sheet"] & (LA < -60)).any() else 0}
    return out, res, geom
