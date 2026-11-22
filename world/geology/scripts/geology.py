"""Layer "geology": rock ages, faults and seismic hazard (experiment E16).

1. Rock age. Every land cell's surface rock (lithology E14) gets a formation age:
   - crystalline rocks (metamorphic, plutonic, shield volcanics) take the age of their province's last tectonic event
     (history E07), a little younger at random (later intrusions and metamorphism);
   - sediments are sampled from the global area-age power law of exposed rocks (area density ~ age^-0.86,
     Wilkinson et al. 2009) between the limits of their setting: platform cover younger than its basement, folded
     sediments older than their orogeny, molasse and rift fill younger than it, basin fill Cenozoic;
   - volcanics of arcs and rifts are as young as their arc or rift; the Doom's pyroclastics ~400 years;
   - surface deposits (alluvium, sand, loess, salt) are Quaternary.
   One smooth random field supplies the quantiles, so that ages form coherent outcrop belts. Ages are grouped into
   the periods of the International Chronostratigraphic Chart.
2. Faults, as one-cell lines with a type and the age of their last movement:
   - plate contacts (plates E06): subduction megathrusts, spreading ridges and transforms, by the boundary types of
     the plate design;
   - contacts of tectonic provinces on continental crust: thrusts at young orogens, arcs and plateaus, normal faults
     at rifts and basins, old thrust fronts and sutures between blocks of different ages. A contact takes the age
     of its younger side's event; active if it moved within 30 Myr, dormant within 66 Myr. Faults of rifted
     margins went dormant when spreading began; continental rifts are active.
3. Seismic hazard (hazard.py, calibrated on Earth): plate boundaries weighted by class and velocity, active
   deforming zones (young orogens, arcs, plateaus, continental rifts) and dormant ones, cut into the four GSHAP
   classes. Rifted margins are passive.

Input:
  --plates, --history, --relief, --ice, --lithology, --canon-relief (npz of those layers)
  world/geology/data/geology.json, world/history/data/history.json, world/plates/data/design.json
Output (in --out):
  geology.npz [6, 1000, 1000]: age_ma float32 (NaN on water and ice), period uint8 (index into periods, 255 none),
                 fault_type uint8 (index into faults.types), fault_age_ma float32 (NaN where no fault),
                 fault_activity uint8 (0 none, 1 active, 2 dormant, 3 inactive),
                 hazard uint8 (index x 100), hazard_class uint8 (0 low, 1 moderate, 2 high, 3 very high)
  stats.json     period shares, the area-age power law against Earth's, fault lengths, hazard by region
  preview_age.png, preview_hazard.png
Deterministic.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
WORLD = HERE.parents[1]
sys.path.insert(0, str(WORLD / "grid" / "scripts"))
sys.path.insert(0, str(WORLD / "climate" / "scripts"))
sys.path.insert(0, str(WORLD / "lithology" / "scripts"))
sys.path.insert(0, str(HERE))
import cubesphere as cs     # noqa: E402
import hazard               # noqa: E402
import lithology as lt      # noqa: E402
import world_climate as wc  # noqa: E402

REGIONS = {"Westeros south of the Wall": (28.0, 63.75, -30.0, 15.0), "Beyond the Wall": (63.75, 90.0, -40.0, 30.0),
           "Western Essos": (25.0, 55.0, 15.0, 50.0)}


def power_law(U, lo, hi, a):
    """Inverse CDF of a density ~ t^-a on [lo, hi] at quantile U."""
    lo = np.maximum(lo, 1e-4)
    hi = np.maximum(hi, lo * 1.0001)
    b = 1 - a
    return (lo ** b + U * (hi ** b - lo ** b)) ** (1 / b)


def rock_ages(rock, deposit, prov, event, fid, order, cfg, U, U2):
    ac = cfg["age"]
    a = ac["power_law_exponent"]
    R = {c: i for i, c in enumerate(order)}
    q, t_min = ac["quaternary_ma"], ac["min_ma"]
    age = np.full(rock.shape, np.nan)
    sed = np.isin(rock, [R["ss"], R["sm"], R["sc"], R["ev"]])
    cryst = np.isin(rock, [R["mt"], R["pa"], R["pi"], R["pb"]])
    young_volc = np.isin(rock, [R["va"], R["vi"], R["py"]])
    vb = rock == R["vb"]
    jitter = 1 - ac["crystalline_jitter"] * U2

    def put(m, v):
        age[m] = np.broadcast_to(v, rock.shape)[m]

    # sediments by setting
    cover = power_law(U, q, event, a)                                  # younger than the basement / the event
    span = ac["folded_sediment_span_ma"]
    folded = power_law(U, event, event + span, a)                      # deposited before the orogeny
    for p in (1, 2, 6):
        put(sed & (prov == p), cover)
    for p in (3, 4):
        post = U2 < ac["post_orogenic_share"][str(p)]
        put(sed & (prov == p) & post, cover)
        put(sed & (prov == p) & ~post, folded)
    put(sed & (prov == 5), power_law(U, q, 100.0, a))                  # fore-arc and back-arc basins
    put(sed & (prov == 7), power_law(U, event, event + ac["plateau_sediment_span_ma"], a))
    put(sed & (prov == 8), power_law(U, t_min, 66.0, a))
    put(sed & np.isnan(age), cover)
    # crystalline rocks: the event that made them
    put(cryst, event * jitter)
    put(cryst & (prov == 5) & (rock != R["mt"]), event + ac["arc_pluton_span_ma"] * U2)
    # volcanics
    put(vb & np.isin(prov, [1, 9]), event * jitter)
    put(vb & np.isin(prov, [2, 3, 7, 8]), cover)                       # flood basalts on platforms
    put(vb & np.isin(prov, [4, 5, 6]), np.maximum(event, q) * U2 + t_min)
    put(young_volc, np.maximum(event, 5.0) * U2 + t_min)
    put(young_volc & (fid == ac["doom_feature"]), ac["doom_age_ma"])
    # unconsolidated sediments and surface deposits
    put(rock == R["su"], power_law(U, t_min, ac["unconsolidated_max_ma"], a))
    put(deposit > 0, power_law(U, t_min, q, a))
    put(np.isin(rock, [R["ig"], R["wb"]]), np.nan)
    return age


def contacts(key, halo):
    """For each cell, a neighbour (4-connected) whose key differs, as (has_contact, neighbour index arrays)."""
    F, J, I = halo
    g = key[F, J, I]
    centre = g[:, 1:-1, 1:-1]
    best = np.zeros(key.shape, bool)
    nb = [(F[:, :-2, 1:-1], J[:, :-2, 1:-1], I[:, :-2, 1:-1]), (F[:, 2:, 1:-1], J[:, 2:, 1:-1], I[:, 2:, 1:-1]),
          (F[:, 1:-1, :-2], J[:, 1:-1, :-2], I[:, 1:-1, :-2]), (F[:, 1:-1, 2:], J[:, 1:-1, 2:], I[:, 1:-1, 2:])]
    pick = [np.zeros(key.shape, int) for _ in range(3)]
    for (f, j, i), gg in zip(nb, (g[:, :-2, 1:-1], g[:, 2:, 1:-1], g[:, 1:-1, :-2], g[:, 1:-1, 2:])):
        m = (gg != centre) & ~best
        for k, arr in enumerate((f, j, i)):
            pick[k][m] = arr[m]
        best |= m
    return best, pick


def main():
    ap = argparse.ArgumentParser()
    for a in ("--plates", "--history", "--relief", "--ice", "--lithology", "--canon-relief", "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((HERE.parent / "data" / "geology.json").read_text())
    hist = json.loads((WORLD / "history" / "data" / "history.json").read_text())
    design = json.loads((WORLD / "plates" / "data" / "design.json").read_text())
    lcfg = json.loads((WORLD / "lithology" / "data" / "lithology.json").read_text())
    n = cs.N
    plate = np.load(args.plates)["plate_cube"].astype(np.int16)
    H = np.load(args.history)
    land = np.load(args.relief)["elevation_m"] > 0
    sheet = np.load(args.ice)["sheet_frac"] == 255
    L = np.load(args.lithology)
    rock, deposit = L["rock"], L["deposit"]
    # provinces and event ages as the lithology drew its rocks from them (through its warp) on land, the history's
    # on the shelves
    prov = np.where(land, L["province"], H["province"]).astype(np.int16)
    event = np.where(land, L["event_age_ma"], H["event_age_ma"]).astype(np.float64)
    crust = H["crust"]
    fid = np.load(args.canon_relief)["feature_id"]
    clat, clon = wc.cell_latlon_deg(n)
    A = np.broadcast_to(cs.cell_area_km2(n), land.shape)
    halo = cs.halo_index(n, 1)

    # ---- 1. rock ages ----
    rng = np.random.default_rng(cfg["age"]["noise"]["seed"])
    U = lt.uniform_field(cfg["age"], n, rng)
    U2 = lt.uniform_field(cfg["age"], n, rng)
    age = rock_ages(rock, deposit, prov, event, fid, lcfg["order"], cfg, U, U2)
    age[~land | sheet] = np.nan
    base = np.array(cfg["periods"]["base_ma"])
    period = np.where(np.isnan(age), 255, np.searchsorted(base, np.nan_to_num(age), side="left")).astype(np.uint8)
    period = np.where(period == len(base), len(base) - 1, period).astype(np.uint8)

    # ---- 2. faults ----
    fc = cfg["faults"]
    T = {t: i for i, t in enumerate(fc["types"])}
    ftype = np.zeros(land.shape, np.uint8)
    fage = np.full(land.shape, np.nan)
    ACT = {a: i for i, a in enumerate(fc["activity"])}
    fact = np.zeros(land.shape, np.uint8)
    # plate contacts: the cell of the lower plate id of each contact pair carries the line
    has, (pf, pj, pi) = contacts(plate, halo)
    other = plate[pf, pj, pi]
    line = has & (plate < other)
    rates = {r["boundary"]: 2 * r["half_rate_km_per_myr"] for r in hist["ridges"]}
    vmed = fc["earth_median_velocity_mm_yr"]
    sources = []
    for b in design["boundaries"]:
        pa, pb = b["between"]
        if pb == 0:
            m = line & ((plate == pa) | (other == pa))
        else:
            m = line & (((plate == pa) & (other == pb)) | ((plate == pb) & (other == pa)))
        if not m.any():
            continue
        cont = crust == 1
        if b["type"] == "trench":
            parts = [("SUB", m, T["subduction megathrust"], vmed["SUB"])]
        elif b["type"] == "ridge":
            v = rates.get(b["name"], vmed["OSR"])
            parts = [("OSR", m & ~cont, T["spreading ridge"], v), ("CRB", m & cont, T["normal"], v)]
        else:
            parts = [("OTF", m & ~cont, T["transform"], vmed["OTF"]), ("CTF", m & cont, T["transform"], vmed["CTF"])]
        for cls, mm, t, v in parts:
            if mm.any():
                ftype[mm] = t
                fage[mm] = 0.0
                fact[mm] = ACT["active"]
                sources.append((cls, v, mm))
    # province contacts on continental crust, drawn on the younger side
    cont = crust == 1
    key = np.where(cont, prov.astype(np.int64) * 100000 + np.round(event * 10).astype(np.int64), -1)
    has, (pf, pj, pi) = contacts(key, halo)
    nprov, nevent, ncont = prov[pf, pj, pi], event[pf, pj, pi], cont[pf, pj, pi]
    edge = has & cont & ncont & (ftype == 0) & ((event < nevent) | ((event == nevent) & (prov > nprov)))
    edge &= ~(np.isin(prov, [6, 8]) & (event < fc["tectonic_min_ma"]))      # post-glacial canon features
    young = prov
    kind = np.select([np.isin(young, [4, 5, 7]), np.isin(young, [6, 8]),
                      (young == 3) & np.isin(nprov, [1, 2])],
                     [T["thrust"], T["normal"], T["thrust"]], T["suture"])
    ftype[edge] = kind[edge]
    fage[edge] = event[edge]
    by_age = np.select([event <= fc["active_max_ma"], event <= fc["dormant_max_ma"]], [ACT["active"], ACT["dormant"]],
                       ACT["inactive"])
    rift_margin = (prov == 6) & (event > fc["continental_rift_max_ma"])
    fact[edge] = np.where(rift_margin, ACT["dormant"], by_age)[edge]

    # ---- 3. seismic hazard ----
    hc = cfg["hazard"]
    reach = 6 * hc["decay_km"]                       # beyond this the index is < 0.3 % of its source's weight
    src = [(cls, v, hazard.distance_km(clat[m], clon[m], clat, clon, reach)) for cls, v, m in sources]
    # deforming zones: young orogens, arcs and plateaus, and continental rifts; rifted margins are passive
    zone_active = cont & ((np.isin(prov, [4, 5, 7]) & (event <= fc["active_max_ma"]))
                          | ((prov == 6) & (event >= fc["tectonic_min_ma"]) & (event <= fc["continental_rift_max_ma"])))
    zone_dormant = cont & np.isin(prov, [4, 5, 7]) & (event > fc["active_max_ma"]) & (event <= fc["dormant_max_ma"])
    zin = np.where(zone_active, hc["zone_weight"], np.where(zone_dormant, fc["dormant_zone_weight"], 0.0))
    zdist = hazard.distance_km(clat[zone_active], clon[zone_active], clat, clon, reach)
    hz = hazard.index(src, zin, zdist, hc)
    hcls = hazard.classes(hz, hc).astype(np.uint8)
    np.savez_compressed(out / "geology.npz", age_ma=age.astype(np.float32), period=period, fault_type=ftype,
                        fault_age_ma=fage.astype(np.float32), fault_activity=fact, hazard=np.round(np.clip(hz, 0, 1) * 100).astype(np.uint8),
                        hazard_class=hcls)

    # ---- statistics ----
    lnd = land & ~sheet & ~np.isnan(age)
    names = cfg["periods"]["names"]
    shares = {nm: round(float(A[lnd & (period == k)].sum() / A[lnd].sum()), 3) for k, nm in enumerate(names)}
    eras = {"Cenozoic": (0, 66.0), "Mesozoic": (66.0, 251.9), "Paleozoic": (251.9, 538.8), "Precambrian": (538.8, 4000)}
    pw = 1 - cfg["age"]["power_law_exponent"]
    era_share = {e: (round(float(A[lnd & (age >= lo) & (age < hi)].sum() / A[lnd].sum()), 3),
                     round((hi / 4000) ** pw - (lo / 4000) ** pw, 3)) for e, (lo, hi) in eras.items()}
    # area-age slope: area per Myr in log bins from 1 Ma to 3000 Ma
    edges = np.logspace(0, np.log10(3000), 15)
    dens = [A[lnd & (age >= e0) & (age < e1)].sum() / (e1 - e0) for e0, e1 in zip(edges[:-1], edges[1:])]
    mid = np.sqrt(edges[:-1] * edges[1:])
    ok = np.array(dens) > 0
    slope = float(np.polyfit(np.log10(mid[ok]), np.log10(np.array(dens)[ok]), 1)[0])
    km = np.sqrt(A)
    flen = {}
    for t, i in T.items():
        if i == 0:
            continue
        m = ftype == i
        flen[t] = {a: round(float(km[m & (fact == k)].sum())) for a, k in ACT.items() if k > 0}
    regions = {}
    for name, (a0, a1, o0, o1) in REGIONS.items():
        m = land & (clat >= a0) & (clat < a1) & (clon >= o0) & (clon < o1)
        t = A[m].sum()
        regions[name] = {"hazard_class_shares": [round(float(A[m & (hcls == k)].sum() / t), 3) for k in range(4)],
                         "period_shares": {nm: v for nm, v in ((nm, round(float(A[m & (period == k)].sum() / t), 3))
                                                               for k, nm in enumerate(names)) if v >= 0.02}}
    import climate as tlayer
    places = {}
    for name, (la, lo) in tlayer.PLACES.items():
        f, j, i = cs.latlon_to_cell(np.radians(la), np.radians(lo), n)
        places[name] = {"rock": lcfg["order"][rock[f, j, i]],
                        "age_ma": None if np.isnan(age[f, j, i]) else round(float(age[f, j, i]), 2),
                        "period": names[period[f, j, i]] if period[f, j, i] < len(names) else None,
                        "hazard_class": int(hcls[f, j, i])}
    stats = {"period_shares_of_land": shares, "era_shares_world_vs_power_law": era_share,
             "area_age_slope": round(slope, 2), "area_age_slope_earth": -0.86, "fault_length_km": flen,
             "hazard_class_shares_land": [round(float(A[land & (hcls == k)].sum() / A[land].sum()), 3) for k in range(4)],
             "regions": regions, "places": places}
    (out / "stats.json").write_text(json.dumps(stats, indent=1, ensure_ascii=False))

    W, Hh = 1440, 720
    plon = -180 + (np.arange(W) + 0.5) * 360 / W
    plat = 90 - (np.arange(Hh) + 0.5) * 180 / Hh
    PLO, PLA = np.meshgrid(plon, plat)
    ff, jj, ii = cs.latlon_to_cell(np.radians(PLA), np.radians(PLO), n)
    pal = np.array(cfg["periods"]["colours"] + [[60, 110, 190]] * (256 - len(names)), np.uint8)
    img = pal[period[ff, jj, ii]]
    img[(ftype[ff, jj, ii] > 0) & land[ff, jj, ii]] = (0, 0, 0)
    Image.fromarray(img).save(out / "preview_age.png")
    hp = np.array([[230, 230, 200], [250, 200, 80], [240, 120, 40], [200, 20, 20]], np.uint8)
    img = np.where(land[ff, jj, ii][..., None], hp[hcls[ff, jj, ii]], np.array([60, 110, 190], np.uint8))
    Image.fromarray(img.astype(np.uint8)).save(out / "preview_hazard.png")
    print(json.dumps(stats, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
