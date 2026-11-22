"""Bioclimatic key to the potential natural vegetation (experiment E15).

The key follows the logic of the BIOME family of models (Prentice et al. 1992; Kaplan et al. 2003): plant types are
admitted by limits of winter cold (mean of the coldest month, Tc), summer heat (growing degree days above 5 and
0 deg C, GDD5 and GDD0) and moisture (alpha = actual / potential evapotranspiration), and the biome is the
vegetation the admitted types make. Wetlands, peat and river floodplains come from the soil and hydrology layers.
Productivity: the Miami model (Lieth 1975). Deterministic.
"""

import numpy as np

DAYS = 365.25 / 12


def bioclim(T, P, PET, whc=150.0, spinup=3):
    """Monthly climate (12, ...) -> dict of Tc, Tw, MAT, GDD5, GDD0, P_year, PET_year, AET_year, alpha, and
    the Mediterranean flag (Koppen Cs: driest month of the warm half-year < 40 mm and < 1/3 of the wettest month
    of the cool half-year), the cool half-year's share of the annual precipitation, and the number of dry months
    (precipitation below half the potential evapotranspiration).
    AET from a Thornthwaite-Mather bucket of `whc` mm (BIOME1 uses 150 mm everywhere), started full and run
    `spinup` years to equilibrium."""
    W = np.full(T.shape[1:], whc, np.float32)
    aet = np.zeros(T.shape[1:], np.float32)
    for yr in range(spinup + 1):
        if yr == spinup:
            aet[:] = 0
        for m in range(12):
            p, e = P[m], PET[m]
            wet = p >= e
            draw = W * (1 - np.exp(-np.maximum(e - p, 0) / whc))
            a = np.where(wet, e, p + draw)
            W = np.where(wet, np.minimum(W + p - e, whc), W - draw)
            if yr == spinup:
                aet += a
    pet_y = PET.sum(0)
    order = np.argsort(T, axis=0)
    warm_half = np.take_along_axis(P, order[6:], 0)
    cool_half = np.take_along_axis(P, order[:6], 0)
    dry = (P < 0.5 * PET).sum(0)
    cool_share = cool_half.sum(0) / np.maximum(P.sum(0), 1.0)
    med = (warm_half.min(0) < 40) & (warm_half.min(0) < cool_half.max(0) / 3)
    return {"Tc": T.min(0), "Tw": T.max(0), "MAT": T.mean(0),
            "GDD5": (np.maximum(T - 5, 0) * DAYS).sum(0), "GDD0": (np.maximum(T, 0) * DAYS).sum(0),
            "P_year": P.sum(0), "PET_year": pet_y, "AET_year": aet,
            "alpha": np.where(pet_y > 1, aet / np.maximum(pet_y, 1), 1.0), "med": med, "cool_share": cool_share,
            "dry_months": dry}


def npp_miami(mat, p_year):
    """Net primary productivity, g dry matter per m2 per year (Lieth 1975): the lesser of the temperature and the
    precipitation limits."""
    nt = 3000 / (1 + np.exp(1.315 - 0.119 * mat))
    npr = 3000 * (1 - np.exp(-0.000664 * p_year))
    return np.minimum(nt, npr)


def classify(inp, cfg):
    """inp: dict of arrays of one shape (climate arrays (12, ...)):
      T (C), P (mm), PET (mm, Thornthwaite), elev_m, lat (deg), land, ice, lake (hydrology code 0-4: 1 lake,
      2 terminal lake, 3 marsh, 4 salt flat), q_m3s (discharge), peat (bool: Histosols of the soil layer).
    Returns (biome codes into cfg['biomes'], the climatic biome before the wetland and river overrides, bioclim)."""
    k = cfg["key"]
    B = {b: i for i, b in enumerate(cfg["biomes"])}
    bc = bioclim(inp["T"], inp["P"], inp["PET"], k["bucket_mm"])
    Tc, Tw, gdd5, gdd0, a = bc["Tc"], bc["Tw"], bc["GDD5"], bc["GDD0"], bc["alpha"]
    elev, lat = inp["elev_m"], np.abs(inp["lat"])
    out = np.full(Tc.shape, -1, np.int16)

    def put(mask, name):
        m = (out < 0) & mask
        out[m] = B[name]

    put(~inp["land"], "Water")
    put(inp["ice"], "Ice")
    put(np.isin(inp["lake"], (1, 2)), "Water")
    alpine = (lat < k["alpine_lat_max"]) & (elev > k["alpine_elev_min_m"])
    # cold limits: no vascular cover, then no trees (tundra), then the open forest-tundra
    put(gdd0 < k["polar_gdd0"], "Polar desert")
    put((gdd5 < k["tree_gdd5"]) & alpine, "Alpine")
    put(gdd5 < k["tree_gdd5"], "Tundra")
    # moisture limits, for every temperature zone
    put(a < k["desert_alpha"], "Desert")
    put(a < k["xeric_alpha"], "Xeric shrubland")
    trop = Tc >= k["tropical_tc"]
    # savanna: moderate moisture, or a long dry season (fire and drought keep the canopy open; Walter's rule)
    put(trop & ((a < k["savanna_alpha"]) | (bc["dry_months"] >= k["savanna_dry_months"])), "Savanna")
    put(trop & (a < k["tropical_moist_alpha"]), "Tropical dry forest")
    put(trop, "Tropical rainforest")
    # extratropics: winter-rain sclerophyll woodland, grassland, then forests by summer heat and winter cold
    winter_rain = bc["med"] & (bc["cool_share"] >= k["med_cool_share"])
    put(winter_rain & (Tc >= k["med_tc"]) & (a >= k["xeric_alpha"]) & (a < k["med_alpha_max"]), "Mediterranean woodland")
    # boreal: short summers, or cold winters (below the broadleaf trees' limit) with moderate summers
    boreal = (gdd5 < k["boreal_gdd5"]) | ((Tc < k["broadleaf_tc_min"]) & (gdd5 < k["boreal_gdd5_cold"]))
    put(boreal & (a >= k["boreal_alpha"]) & (gdd5 < k["forest_tundra_gdd5"]), "Forest-tundra")
    put(boreal & (a >= k["boreal_alpha"]), "Boreal forest")
    put(a < k["forest_alpha"], "Temperate grassland")
    put(boreal, "Temperate grassland")
    put((Tc >= k["warm_temperate_tc"]) & (gdd5 >= k["warm_temperate_gdd5"]), "Subtropical evergreen forest")
    # conifers: montane forests, and the temperate rainforests of cool, very wet oceanic coasts
    rainforest = (bc["P_year"] >= k["temperate_rainforest_p_mm"]) & (Tw < k["temperate_rainforest_tw_max"])
    put((elev > k["conifer_elev_m"]) | rainforest, "Temperate conifer forest")
    put((gdd5 < k["mixed_gdd5"]) | (Tc < k["broadleaf_tc_min"]), "Temperate mixed forest")
    put(np.ones(Tc.shape, bool), "Temperate broadleaf forest")
    # high plateaus and mountains above the forest: montane grassland and alpine meadow
    high = (lat < k["alpine_lat_max"]) & (elev > k["montane_elev_m"])
    open_land = np.isin(out, [B["Tundra"], B["Temperate grassland"], B["Xeric shrubland"], B["Savanna"]])
    out = np.where(high & open_land, B["Alpine"], out)
    climatic = out.copy()
    # water and soils: marshes and peat bogs (swamp forest where the climate carries forest), river galleries
    # through dry land
    forest = np.isin(climatic, [B[b] for b in cfg["forest_biomes"]])
    wet = inp["land"] & ~inp["ice"] & (climatic != B["Water"]) & (climatic != B["Polar desert"]) \
        & ((inp["lake"] == 3) | inp["peat"])
    swamp = wet & forest & (gdd5 >= k["swamp_forest_gdd5"])
    out = np.where(swamp, B["Swamp forest"], np.where(wet, B["Wetland"], out))
    dry = np.isin(climatic, [B["Desert"], B["Xeric shrubland"], B["Temperate grassland"], B["Savanna"]])
    out = np.where(dry & (inp["q_m3s"] >= k["riparian_q_m3s"]) & ~wet, B["Riparian woodland"], out)
    out = np.where(inp["lake"] == 4, B["Salt flat"], out)
    return out, climatic, bc
