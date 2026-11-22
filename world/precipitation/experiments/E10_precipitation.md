# E10 — Precipitation: a moisture-budget model, tested on Earth

Status: **done** — run m is the layer. The Dorne conflict (runs i–j) was resolved by adding the monsoon (owner decision 2).

## Goal

Monthly precipitation for every cell, from physics, built on the temperature layer (E09) and the relief (E08), and
tested on Earth before it is run on the world. Canon checks: the named deserts (Dornish desert, Great Sand Sea, Red
Waste, Grey Waste) should be arid (< ~250 mm/year); Westeros as a whole is temperate and green.

## Model (scripts/moisture.py)

A steady, vertically integrated water-vapour budget on a latitude–longitude grid, solved each month:

    div(W V) − div(K ∇W) + k W = E,        P = k W

- **W** column water vapour (mm), **V** low-level wind, **K** eddy diffusivity (transient eddies carry moisture down
  the gradient; Pierrehumbert 2002, Frierson et al. 2006 use K ~ 10⁶ m²/s), **k** precipitation rate (1/s),
  **E** evaporation. P = E globally by construction (mass conservation).
- **Winds** (`winds`): zonal belts (trades, westerlies stronger in winter, polar easterlies); Hadley / Ferrel /
  polar meridional cells about an ITCZ that moves with the season and, optionally, towards the local thermal
  equator; a thermal component towards warm ground and its geostrophic counterpart (anticlockwise round warm lows
  in the north: monsoon inflow, subtropical highs over cool seas).
- **Precipitation rate k:** a base time scale τ0, longer in warm air (Clausius–Clapeyron), multiplied by
  exp(boost), with boost the sum of:
  - low-level convergence (ascent);
  - orographic upslope wind;
  - **Hadley subsidence:** −s·G(|x|), a Gaussian band about 23° from the local ITCZ, stronger on land;
  - **Sverdrup ascent/descent:** β v = f ∂w/∂z in the lower troposphere, so poleward flow rises and equatorward flow
    sinks (Rodwell & Hoskins 2001): the wet western and dry eastern flanks of the subtropical highs.
- **Monsoon (Gill circulation):** the damped linear shallow-water response of the first baroclinic mode (Gill 1980,
  on the sphere: `gill_operator`) to a heating ∝ the month's temperature anomaly from the zonal mean, faded
  poleward of 35°. It gives low-level inflow into the summer continents, westerlies on the west side of the heating
  and the cross-equatorial jet, and the matching outflow from cool winter land. Its wind is added to the analytic
  winds.
- **Evaporation:** over open sea, bulk formula ∝ (wind + 1)·q_sat(T); none over sea ice; over land the smaller of a
  potential rate ∝ q_sat(T_ground) (zero below 0 °C) and a fraction of last pass's local precipitation.
- **Numerics:** first-order upwind finite volumes, sparse direct solve per month, two passes for land evaporation.
  The model is calibrated and run at **2°**: upwind numerical diffusion depends on the grid step, so a 1° run is
  not the calibrated model (tested: at 1° the world's land is markedly drier).

## Earth test (scripts/earth_precip_validation.py)

- Inputs: Earth temperature and sea ice from the E09 model itself (so the test sees the same temperature errors
  as the world), Earth land fraction (global-land-mask), ETOPO1 10′ topography (fatiando-data release,
  DOI 10.5281/zenodo.5882203), averaged to the grid.
- 109 stations (`data/earth_precip_stations.json`, from Wikipedia climate tables: annual, January, July), sampled
  from the land cells around each station.
- Approximate GPCP zonal means (`data/earth_zonal_gpcp.json`, Adler et al. 2003), 15 latitudes.
- Scores: RMSE of log10(model/station) annual; share within a factor 2; seasonal contrast (Jan − Jul)/year; bias of
  the arid stations (< 250 mm); zonal-mean log RMSE; global mean (Earth ~990 mm).

## Runs

| Run | Change | log RMSE | within ×2 | arid bias (log) | zonal log RMSE | Notes |
|---|---|---|---|---|---|---|
| a | first version, 1°, Nelder–Mead calibration | 0.43 | 59 % | ≈ +0.5 | — | deserts ~5× too wet (Cairo, Muscat); Dornish desert 759 mm |
| b | precipitation ∝ exp(a·(W/W_sat − r0)) (Bretherton et al. 2004), fixed-point | — | — | — | — | first form broke conservation (P from the damped W): global 49 000 mm. Fixed by a final solve |
| c | b, conserving, calibrated at 2° | 0.44 | 57 % | ≈ +0.45 | — | calibration drove a → 0: the humidity dependence did not help; switched off (`rh_iterations` 0) |
| d | stations sampled from land cells only | 0.44 | 56 % | — | — | no change: the coastal-desert error is in the land cells themselves |
| e | ITCZ towards the local thermal equator (Schneider et al. 2014) | 0.45 | 54 % | — | — | little effect; over the Sahara the thermal equator is the heat low, north of the real rain belt |
| f | + Hadley subsidence band | 0.52 | 43 % | +0.30 | — | Cairo 61 mm, Tamanrasset 62 mm; but Miami, Delhi far too dry |
| g | + Sverdrup term, calibrated (18 parameters) | 0.45 | 49 % | +0.22 | 0.30 | the fit drove K to 2·10⁴ m²/s and the westerlies to peak in summer: **no mid-latitude storm-track rain** (50° N ocean ~430 mm vs ~1000). World: Westeros 48 % below 250 mm, Winterfell 200 mm, King's Landing 280 mm — rejected |
| h | K fixed at 10⁶ m²/s (literature) | 0.45 | 55 % | — | — | better on stations too; the zonal profile gains its storm tracks |
| i | + zonal-mean GPCP term in the objective, recalibrated (2 starts) | 0.45 | 56 % | +0.36 | 0.08 | Dornish desert 599 mm: canon conflict, see below |
| j | statistical station correction | 0.43 | — | +0.6 | — | failed, see below |
| k | monsoon as a local Ekman-balanced heat-low wind | 0.47–0.79 | 30–55 % | +0.27–0.37 | — | a Somali jet appears but turns north into Arabia; worse |
| l | monsoon as a Gill circulation driven by last pass's rain | 0.73–1.05 | 26–34 % | ≈ 0 | — | runaway: the heating sits over the ocean ITCZ and pulls the rain off the land |
| m0 | Gill circulation driven by the temperature anomaly, calibrated | 0.39 | 54 % | +0.16 | 0.08 | Lisbon 3800, Seattle 2200 mm: the forcing by cold winter continents makes a spurious mid-latitude south-westerly |
| **m** | m0 with the Gill forcing faded poleward of 35°, recalibrated | **0.39** | **56 %** | **+0.22** | **0.09** | chosen (`data/precipitation.json`) |

Calibration notes for run m (Rule 1 — where parameters are effective rather than literal):
- `hadley_ms` → 0: the mean meridional moisture transport is carried by the eddy diffusion (K ≈ 3·10⁶ m²/s, upper
  end of the literature range) rather than by the analytic Hadley wind.
- `land_recycling` ≈ 0.11: an effective parameter; with no soil-moisture store it is not the observed land ET/P
  (~0.6).
- `westerly_winter_boost` ≈ 0.06: small but of the physical sign (it was bounded at ≥ 0 after run g).
- `boost_max` ≈ 7.5: the cap on the rain-rate multiplier is rarely reached; a cap of 3–4 was tested and is worse.
- Calibrations ran in the background with a resumable Nelder–Mead (every evaluation logged; a restart builds its
  simplex from the best points), because the sandbox stops background processes.

### Earth skill of run m

- 109 stations: log RMSE 0.39 (typical error ×2.5), 56 % within a factor 2, seasonal-contrast RMSE 0.12.
- Zonal means within ~20 % of GPCP: 367, 741, 905, 695, 578, 1404, 1888, 2204, 1200, 898, 751, 1084, 992, 516,
  276 mm at 70° N … 70° S; global mean 1026 mm.
- Arid stations (14, < 250 mm): mean bias +0.22 log (×1.7 too wet). Peninsular deserts improved from run i:
  Doha 234 / 75 (was 592), Muscat 220 / 100 (was 532); Cairo 140 / 25, Tamanrasset 48 / 58.
- Mid-latitudes: London 770 / 602, Madrid 341 / 423, Rome 555 / 799, Moscow 451 / 714, Chicago 318 / 1038.
- **Known failures** (model / station, mm):
  - monsoons still too dry: Mumbai 449 / 2213, Delhi 66 / 759, Sylhet 894 / 7377 (the onset is right in season,
    July wet and January dry, but the amount is 3–8× low);
  - some windward west coasts too wet in winter: Seattle 2822 / 945, Lisbon 1551 / 774;
  - the eastern United States too dry (Chicago, Miami ×3).

## World (run m, 2°)

| Place | Year | January | July |
|---|---|---|---|
| Castle Black | 701 | 50 | 72 |
| Winterfell | 690 | 66 | 39 |
| King's Landing | 286 | 26 | 16 |
| Oldtown | 340 | 85 | 3 |
| Sunspear | 221 | 36 | 7 |
| Volantis | 479 | 154 | 1 |
| Qarth | 382 | 59 | 11 |
| Asshai | 1317 | 35 | 126 |

Global mean 985 mm; land 778 mm; Westeros below 250 mm: 17 %. The south of Westeros and the Free Cities have a
Mediterranean regime (wet winter, dry summer).

| Canon desert | Median | 10th–90th percentile | Share < 250 mm | Share < 400 mm |
|---|---|---|---|---|
| Red Waste | 85 | 33–163 | 98 % | 100 % |
| Grey Waste and Cannibal Sands | 133 | 87–243 | 91 % | 100 % |
| Great Sand Sea | 169 | 94–258 | 87 % | 100 % |
| Dornish desert | 268 | 175–374 | 41 % | 95 % |

All four canon deserts are arid or semi-arid. Allowing for the model's Earth bias over arid land (×1.7 too wet), the
Dornish desert is likely below 200 mm. **Canon check passed.** Remark: King's Landing (286 mm) is dry for the
Crownlands; Earth analogue Madrid is modelled 341 / 423 mm, so the model is not obviously wet-biased there. Not a
canon conflict (canon gives no rainfall), noted for the vegetation layer.

## Addendum (2026-11-22, from E11): snow on high ground

The temperature model now puts snow where the ground, not sea level, is below freezing (E09 addendum). Earth's
climate became cooler at high latitudes and the global precipitation fell from 1026 to 944 mm. `ocean_e0_mm_day` was
refitted alone to the global mean, 4.55 → 4.77. Earth test after the refit:
- log RMSE 0.39 and 58 % within a factor 2 (both unchanged);
- arid bias +0.25 log;
- zonal log RMSE 0.11;
- global mean 989 mm.

The world is now run through `world/climate/scripts/world_climate.py` and takes an optional `--ice`.

## Addendum (2026-11-22, from E12): recalibration with Earth rivers, reverted

The model was refitted with Earth's river discharge and runoff to the sea added to the objective (E12, run p10).
Earth improved a little: stations 0.388 → 0.376, rivers 0.71 → 0.63. The Reach stayed as dry as before. The wetter
North grew a northern ice sheet in the ice layer again. Reverted; details in E12.

## Canon conflict (Rule 3), resolved

The Dornish desert comes out at ~600 mm — semi-humid to semi-arid, not a desert. Dorne is a peninsula in a warm sea,
which is exactly the setting where the model fails on Earth (Doha, Muscat ×5–8 too wet).

**Owner decision 1 (2026-11-21): try an Earth-fitted statistical correction**, as for temperature.

### Run j — station correction (failed)

log10(station/model) regressed on predictors the world run also has, every combination of 1–3 of: annual
subsidence index, sea fraction around the cell, their product, Sverdrup vertical velocity, model log P, seasonal
temperature amplitude, annual and warmest-month temperature. Leave-one-out over 109 stations:

| Predictors | LOO log RMSE | arid bias after |
|---|---|---|
| none | 0.446 | +0.36 |
| model log P | 0.432 | +0.64 |
| subsidence | 0.445 | +0.52 |
| sea fraction, log P | 0.432 | +0.63 |

No combination improves the fit by more than 0.014, and every one makes the arid stations worse. The reason: the two
large errors sit in the same regime with opposite signs. Monsoon stations (Mumbai, Delhi, Miami: model far too dry)
and peninsular deserts (Doha, Muscat: far too wet) have the same subsidence, coast and temperature predictors. A
regression on them cancels out. The likely common cause is physical: the model has no strong summer monsoon
inflow, so the moisture of the Arabian Sea that on Earth feeds India stays and rains over Arabia.

Back to the owner (Rule 3). **Owner decision 2 (2026-11-21): add the monsoon** — runs k–m above.
