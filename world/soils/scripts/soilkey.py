"""Climate and terrain key to the WRB Reference Soil Groups (experiment E13).

Soils are named by the factor that most clearly conditions them, in the order of the sets of the FAO lecture notes
(World Soil Resources Report 94, 2001): organic and waterlogged soils, then soils conditioned by topography, parent
material and permafrost, then by climate (arid, steppic, wet tropical, temperate and boreal). Every rule uses
inputs that the world's layers provide: monthly temperature and precipitation, potential evapotranspiration, local
relief, height, and the hydrology layer's lakes, marshes, salt flats and discharge. Parent material is known only
where it is volcanic (Andosols). Deterministic.
"""

import numpy as np


def dry_months(P_month, PET_month):
    """Months whose precipitation is below half the potential evapotranspiration (a dry season month)."""
    return (P_month < 0.5 * PET_month).sum(0)


def moisture_class(P_year, PET_year):
    """Thornthwaite (1948) moisture index Im = 100 (P / PET - 1) and its classes:
    0 arid (Im < -40), 1 semi-arid (-40..-20), 2 dry subhumid (-20..0), 3 moist subhumid (0..20),
    4 humid (20..100), 5 perhumid (> 100)."""
    im = 100 * (P_year / np.maximum(PET_year, 1.0) - 1)
    return np.digitize(im, [-40, -20, 0, 20, 100]), im


def classify(inp, cfg):
    """inp: dict of arrays of one shape:
      T (12, ...) C, P (12, ...) mm, PET (12, ...) mm (Thornthwaite), relief_m (local relief over ~50 km), elev_m,
      land, ice, lake (hydrology code 0-4), q_m3s (discharge), volcanic (bool), basin (bool: closed depression).
    Returns int codes into cfg['groups']."""
    k = cfg["key"]
    G = {g: i for i, g in enumerate(cfg["groups"])}
    T, P, PET = inp["T"], inp["P"], inp["PET"]
    mat = T.mean(0)
    tcold = T.min(0)
    twarm = T.max(0)
    Pa = P.sum(0)
    mc, im = moisture_class(Pa, PET.sum(0))
    ARID, SEMI, DRYSUB, MOIST, HUMID, PERHUMID = range(6)
    dry = dry_months(P, PET)
    rel = inp["relief_m"]
    elev = inp["elev_m"]
    lake = inp["lake"]
    basin = inp.get("basin", np.zeros(mat.shape, bool))
    lowland = basin | ((elev < k["lowland_elev_m"]) & (rel < k["lowland_relief_m"]))
    # lithology (E14), when given: rock classes of GLiM and surface deposits
    rock = inp.get("rock")
    if rock is not None:
        R = {c: i for i, c in enumerate(cfg["rock_order"])}
        dep = inp["deposit"]
        young_volc = np.isin(rock, [R["vi"], R["va"], R["py"]])
        basic = np.isin(inp["bedrock"], [R["vb"], R["pb"]])
        carbonate = inp["bedrock"] == R["sc"]
        sand = dep == 3
        loess = dep == 4
        # smectite-forming: basalts, alluvial and basin clays, shales and marls of mixed sediments on low ground
        clayey = basic | (dep == 1) | basin | ((inp["bedrock"] == R["sm"]) & lowland)
    else:
        young_volc, basic, carbonate, sand, loess = (np.zeros(mat.shape, bool),) * 5
        clayey = lowland
    out = np.full(mat.shape, -1, np.int16)

    def put(mask, name):
        m = (out < 0) & mask
        out[m] = G[name]

    put(~inp["land"], "Water")
    put(inp["ice"], "Ice")
    put(np.isin(lake, (1, 2)), "Water")
    marsh = lake == 3
    cool = mat < k["boreal_mat_C"]
    # organic soils: marshes in cool climates, and wet flats of the boreal and arctic lowlands (peat also forms over
    # permafrost, as in the West Siberian lowland)
    put(marsh & (mat < k["histosol_mat_max_C"]) & (mc >= MOIST), "Histosols")
    put(cool & (mc >= HUMID) & (rel < k["boreal_wet_flat_relief_m"]), "Histosols")
    put(lake == 4, "Solonchaks")
    put(mat <= k["permafrost_mat_C"], "Cryosols")
    # topography and parent material
    put(rel > k["leptosol_relief_m"], "Leptosols")
    put((mc == ARID) & (rel > k["leptosol_arid_relief_m"]), "Leptosols")
    # ash soils need young volcanic ejecta and a moist climate (FAO: humid volcanic regions)
    put(inp["volcanic"] | (young_volc & (mc >= MOIST)), "Andosols")
    put(inp["q_m3s"] >= k["fluvisol_q_m3s"], "Fluvisols")
    put(marsh, "Gleysols")
    # deserts (FAO: < ~250 mm; Thornthwaite arid)
    desert = (mc == ARID) & (Pa < k["desert_p_mm"])
    put(desert & basin & (Pa < k["gypsisol_p_mm"]), "Gypsisols")
    put(desert & sand, "Arenosols")
    put(desert & (rel < k["arenosol_relief_max_m"]) & (tcold >= k["arenosol_cold_month_C"]) & (rock is None),
        "Arenosols")
    put(desert & (Pa < k["hyperarid_p_mm"]), "Leptosols")
    put(desert, "Calcisols")
    # old, stable surfaces of the humid tropics
    trop = tcold >= k["tropical_cold_month_C"]
    # basic rocks of the humid tropics: deep, well-structured red clays
    put(trop & (mc >= MOIST) & basic & (rel <= k["ferralsol_relief_max_m"] * 1.5), "Nitisols")
    put(trop & (mc >= MOIST) & (rel <= k["ferralsol_relief_max_m"]) & (dry <= k["dry_months_ferralsol_max"]),
        "Ferralsols")
    # seasonal clay plains: smectitic sediments of basins and lowlands, wet and dry seasons
    warm = tcold >= k["warm_cold_month_C"]
    put(warm & clayey & (rel < k["lowland_relief_m"] * 2) & (dry >= k["vertisol_dry_months_min"]) &
        np.isin(mc, (SEMI, DRYSUB, MOIST)), "Vertisols")
    # sandy plains of the dry tropics (Sahel, Kalahari): summer rains (the wettest month in the warm half year);
    # dry plains with winter rains are Mediterranean (Luvisols, Calcisols)
    summer_rain = np.take_along_axis(T, P.argmax(0)[None], 0)[0] >= mat
    put(warm & (mc <= SEMI) & sand, "Arenosols")
    put(warm & (mc <= SEMI) & (rel < k["arenosol_relief_max_m"]) & (tcold >= k["arenosol_cold_month_C"]) &
        summer_rain & (rock is None), "Arenosols")
    # limestone hills of humid lands: thin rendzinas
    put(carbonate & (mc >= MOIST) & (rel >= k["cambisol_relief_m"]), "Leptosols")
    # steppes: cold winters (FAO: Kastanozems 200-400 mm, Chernozems ~450-500 mm with P ~ evaporation,
    # Phaeozems where there is a surplus)
    cold = (tcold < k["steppe_cold_month_C"]) & (mat > k["steppe_mat_min_C"])
    put(cold & (mc <= SEMI), "Kastanozems")
    put(cold & (mc == DRYSUB), "Chernozems")
    put(cold & (mc == MOIST) & (twarm >= k["phaeozem_warm_month_C"]), "Phaeozems")
    # wet tropics
    put(trop & (elev > k["highland_tropics_m"]) & (mc >= HUMID), "Umbrisols")
    put(trop & (mc >= MOIST), "Acrisols")
    put(trop, "Lixisols")
    # subtropics, warm temperate and Mediterranean
    put(warm & (mc <= SEMI), "Calcisols")
    put(warm & (mc >= HUMID), "Acrisols")
    put(warm & (rel >= k["cambisol_relief_m"]), "Cambisols")
    put(warm, "Luvisols")
    # cool temperate and boreal
    put(cool & (twarm < k["boreal_warm_month_C"]) & (mc >= MOIST), "Podzols")
    put(sand & (mc >= MOIST), "Podzols")             # leached quartz sands of humid lands
    put(loess & (mc >= MOIST), "Luvisols")           # the loess belt of humid temperate lands
    put((tcold < k["albeluvisol_cold_month_C"]) & (mat < k["albeluvisol_mat_max_C"]) & (mc >= MOIST),
        "Albeluvisols")
    put((mc == PERHUMID) & (tcold >= k["umbrisol_cold_month_C"]), "Umbrisols")
    put(rel >= k["cambisol_relief_m"], "Cambisols")
    put(mc <= SEMI, "Kastanozems")
    put(mc <= HUMID, "Luvisols")
    put(np.ones(mat.shape, bool), "Cambisols")
    return out
