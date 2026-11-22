# E09 — Temperature: a seasonal energy balance model, tested on Earth

## Goal

Monthly mean near-surface temperature for every cell, from physics the owner chose (2026-11-21): a seasonal energy
balance model on the sphere, calibrated and tested on Earth before it is run on the world. The first canon check is
whether the Haunted Forest (64–66° N) has a summer warm enough for trees: the polar tree line follows a warmest
month of about 10 °C.

## Model (scripts/ebm.py)

- **Grid and state.** One temperature per cell of a 2° latitude–longitude grid, for the surface and lower air
  together, reduced to sea level.
- **Sunlight.** Daily insolation from Earth's orbit.
- **Albedo.** Planetary albedo a0 + a2·P2(sin φ); ice and snow albedo 0.62.
- **Outgoing radiation.** Linear outgoing longwave radiation A + B·T (Budyko/North).
- **Heat transport.** Diffusion on the sphere (North 1975), plus advection of the air's heat capacity by a zonal
  wind profile (trades, westerlies, polar easterlies).
- **Heat capacity.** A 20–25 m ocean mixed layer, and a small land heat capacity.
- **Sea ice.** Zero-layer thermodynamic sea ice (Semtner 1976): it holds the surface at 0 °C while it melts, and grows
  by conduction.
- **Continentality.** An index M from the distance to the sea, measured upwind along the wind. Far from the sea
  there are fewer clouds (lower albedo) and drier air (higher outgoing radiation).
- **Time stepping.** Backward Euler, one day per step, sparse LU solves; 12 years of spin-up.

Tried and switched off: wind-driven ocean gyres (Stommel 1948 streamfunction on the actual coasts, heat carried by
the depth-integrated transport). The gyres came out right in shape, with western boundary currents, but at ~4 Sv
instead of ~30 Sv. Scaled up, they did not warm north-west Europe and they over-warmed Vladivostok. The code is kept
with `heat_factor` 0.

## Earth test (scripts/earth_validation.py)

The test uses Earth's land mask (pip package global-land-mask) and 78 stations from "List of cities by average
temperature" (Wikipedia), with January, July and annual means. Station values are reduced to sea level at
6.5 K/km. The free parameters, 11 of them, were fitted by Nelder–Mead:
- radiation constants and albedo;
- diffusion and zonal diffusion;
- westerly strength;
- land heat capacity and ocean mixed-layer depth;
- the three continentality parameters.

Where the fit wandered to unphysical values (land heat capacity at its bound, continentality scale ~100 km), an
earlier, physically plausible point was kept.

**Owner decision (CONSTRAINTS decision 7).** The model has no thermohaline overturning. Europe north of 45° N
(17 stations) is ~12 °C too cold in winter, because there the North Atlantic Current does work the model cannot. The
owner decided that the Sunset Sea has no such current. Those stations are excluded from calibration and reported
separately; the Pacific coasts are the analogue.

| Stage | Stations | RMSE Jan / Jul / year, °C |
|---|---|---|
| first run, textbook parameters | 68 | 8.7 / 6.4 / 4.1 |
| calibrated | 68 | 5.1 / 4.2 / 3.5 |
| + polar and highland stations | 78 | 5.6 / 4.9 / 4.1 |
| + sea ice, European stations excluded, recalibrated | 61 | **4.9 / 4.5 / 3.5** |
| + statistical correction (below), leave-one-out | 61 | **4.5 (winter) / 4.3 (summer)** |

**Statistical correction.** The residual (model − station) is fitted by least squares for the winter month and the
summer month: b0 + b1·z + b2·M, with z the elevation in km (capped at 2.5 km) and M the continentality index. On land
the fitted residual is subtracted, with the coefficients moving between winter and summer through the year. It
captures the warmth of high plateaus (an effective lapse rate of ~3.5 K/km rather than 6.5) and some of the winter
cold of continental interiors.

## Known errors on Earth (after calibration)

- **Tropics too warm by ~3 °C.** The model has no evaporation or deep convection; the equatorial zonal mean is 31 °C
  against ~27 °C.
- **Global mean 16.8 °C** against 14 °C (NASA, 1951–1980).
- **Continental interiors swing too little.** Winters are too warm (Yakutsk +15, Winnipeg +9.5, Yellowknife
  +8.4 °C) and summers too cold (Fairbanks −7.6, Yellowknife −6.2, Anchorage −3.9 °C). Land is coupled too strongly
  to the ocean by horizontal mixing.
- **Pacific north-west too cold in winter by ~5 °C** (Vancouver, Seattle, Anchorage).
- **Hot deserts and plateaus too cold in summer by 10–13 °C** (Tehran, Riyadh, Denver), partly corrected.

## First run on the world (run w1)

| Place | Elevation, m | Jan / Jul / year, °C |
|---|---|---|
| Castle Black | 471 | −17.6 / 4.9 / −5.9 |
| Winterfell | 295 | −11.5 / 8.7 / −1.4 |
| King's Landing | 297 | 4.6 / 18.7 / 11.6 |
| Oldtown | 0 | 15.7 / 23.9 / 19.9 |
| Sunspear | 24 | 16.3 / 25.0 / 20.8 |
| Volantis | 246 | 14.3 / 26.0 / 20.4 |
| Qarth | 76 | 17.6 / 31.2 / 25.0 |

- **Haunted Forest:** warmest month 5.5 °C on average, so no cell reaches 10 °C. That would make it tundra, not
  forest.
- **Land of Always Winter:** warmest month −0.1 °C; 53 % stays below 0 °C all year.

**Verdict: not accepted.** The north of Westeros falls in exactly the regime where the model is weakest on Earth:
summers in subarctic land are 4–8 °C too cold (Fairbanks, Yellowknife, Anchorage). A Haunted Forest July of 5.5 °C
plus that bias gives ~10–13 °C, which is forest. So the canon check cannot be decided with this model, and the
coupling of land and air has to be improved first.

## Two-layer model (runs w2–w5)

`run_two_layer` in ebm.py keeps the surface (land, mixed-layer ocean or sea ice) and the air as two temperatures.
- Only the air is mixed and carried by the wind.
- The surface exchanges heat with the air: strongly when the surface is warmer (convection), weakly when it is
  colder (stable inversions).
- The surface also loses heat through the atmospheric window.
- Sea ice became a zero-layer ice surface over water held at freezing point. A first version changed the heat
  capacity of melting cells only once a month, which overheated fresh melt water (40 °C stripes at the ice edge)
  and gave the calibration false heat; it was replaced.
- Calibrated again (14 parameters; an earlier point with a plausible land heat capacity of 2.2 × 10⁶ J m⁻² K⁻¹
  was kept).

| | One layer | Two layers |
|---|---|---|
| RMSE Jan / Jul / year (61 stations) | 4.9 / 4.5 / 3.5 | 5.3 / 4.2 / 3.7 |
| After correction, leave-one-out winter / summer | 4.5 / 4.3 | 4.7 / 3.6 |
| Fairbanks July | −7.6 | −2.6 |
| Yellowknife July | −6.2 | −1.1 |
| Anchorage July | −3.9 | +0.5 |
| Arkhangelsk July | +0.5 | −1.3 |

The two-layer model fixes subarctic summers on land, the regime that matters for the North of Westeros. Interior
winters stay too warm (Yakutsk +18 °C).

**Sea ice is the weak point.** NSIDC's 1981–2010 averages are Arctic 15.6 M km² in March and ~6.5 in September,
Antarctic 18.6 in September and ~2.5 at summer's end. The model gives the Arctic 26 / 0 and the Antarctic 41 / 3.5:
too much winter ice and no perennial Arctic ice. A dedicated sea-ice fit (ice albedo, ocean heat under the ice)
traded station accuracy for extent and was rejected.

## World with the two-layer model (run w5)

| Place | Jan / Jul / year, °C |
|---|---|
| Castle Black | −21.6 / 6.5 / −7.1 |
| Winterfell | −14.0 / 9.4 / −2.5 |
| King's Landing | 1.8 / 18.4 / 10.3 |
| Oldtown | 15.8 / 24.0 / 20.2 |
| Sunspear | 14.0 / 24.6 / 19.7 |
| Volantis | 12.7 / 25.5 / 19.6 |
| Qarth | 17.2 / 29.5 / 23.9 |

- **Haunted Forest:** warmest month 6.1 °C; no cell reaches 10 °C.
- **Land of Always Winter:** warmest month 2.6 °C.
- **Sea ice:** the north polar sea keeps 22 M km² of ice in September (Earth's Arctic ~6.5). Its July edge reaches
  65° N on both sides of Westeros. The world's north polar region is mostly ocean, while Earth's Arctic Ocean is
  ringed by land that warms in summer. So a cold, icy northern sea is a real consequence of the geography. Because
  the model melts sea ice too easily on Earth, the true extent may be larger still.

**Verdict.** Without a warm current (decision 7), the North of Westeros is like Labrador and Baffin Island:
- Winterfell is at about the tree line;
- the Haunted Forest is forest-tundra rather than forest;
- the Land of Always Winter is above freezing in summer at sea level, and would be "always winter" only under an
  ice sheet (a later layer).

This contradicts canon (WE7, the Haunted Forest, the Wolfswood). Under rule 3 the work stops here for the owner's
decision.

## Options for the owner

1. **Keep decision 7** and log a canon deviation: the Haunted Forest is sparse forest-tundra, and the North is
   forested only in its south.
2. **Reverse decision 7: a warm drift in the north-east Sunset Sea.** The heat it brings would be a prescribed
   ocean heat flux (the "Q-flux" of slab-ocean models), calibrated on the 17 European stations excluded so far. The
   North becomes Norway-like, with forest to the Wall and beyond.
3. **Move the map south** (change decision 3). Every layer is regenerated by script, but the canon latitudes of the
   climate bands change.

## Owner's choice: the warm drift (run w6)

The owner chose option 2 (2026-11-21), reversing decision 7. The warm drift is a Q-flux, a heat flux added to the
mixed layer of one basin's ocean:
- Q0·exp(−((φ − φ0)/w)²)·exp(−d/L), where d is the distance to the coast to the east along the parallel;
- the same heat is taken from the basin's ocean between 0 and 30° N, so the field sums to zero.

Calibrated on the North Atlantic with all 78 stations (the other parameters fixed): Q0 = 88 W m⁻², φ0 = 54° N,
w = 15°, L = 1860 km. The RMSE is 5.4 / 4.0 / 3.6 °C. After correction, leave-one-out gives 5.0 (winter) and
3.2 °C (summer). The European stations improve most in winter: Bergen goes from −13 to −6 °C, Reykjavik from
−12 to −6 °C.

On the world the same field lies in the northern Sunset Sea (`world_qflux_box`, 170° E to 15° W).

| Place | Jan / Jul / year, °C |
|---|---|
| Castle Black | −18.2 / 8.2 / −3.8 |
| Winterfell | −9.0 / 11.3 / 0.4 |
| King's Landing | 3.8 / 19.6 / 11.9 |
| Oldtown | 17.1 / 25.2 / 21.4 |
| Sunspear | 15.8 / 25.2 / 20.9 |
| Volantis | 14.4 / 26.0 / 20.7 |
| Qarth | 18.6 / 29.7 / 24.8 |

- **Haunted Forest:** warmest month 7.6 °C (maximum 8.6 °C).
- **Land of Always Winter:** warmest month 3.1 °C.
- **Sea ice:** northern sea ice falls to 15 M km² in July, with its edge at 73–83° N near Westeros.

The drift warms the North's winters a lot but its summers by only 1.5–2 °C, as on Earth (Bergen's July is 15.5 °C,
Tromsø's 11.8 °C). Land beyond the Wall at 64–70° N has July means of 5–9 °C:
- comparable to Iceland, rather than to Finland;
- the land at 64° N is ~1800 km wide but lies between two oceans with ice-covered polar sea and the Land of Always
  Winter to the north;
- on Earth the model gets Scandinavian summers within 1–2 °C (Helsinki −0.8 °C, Arkhangelsk −1.3 °C), so this is a
  consequence of the world's geography rather than a model error.

Uncertainty for the North's summers is about ±3 °C (the model's subarctic errors on Earth). The Haunted Forest's
7.6 °C is below the ~10 °C tree line even at the upper end, so it is forest-tundra or sparse woodland in sheltered
valleys, as in Iceland or the Kola peninsula. Winterfell (11.3 °C) is just above the tree line.

**Decision (owner, 2026-11-21).** Run w6 is accepted as the temperature layer. The Haunted Forest as sparse taiga at
the edge of the forest-tundra is logged as canon deviation 5.

Still open, for later layers:
- ice sheets, which would make the Land of Always Winter "always winter";
- precipitation;
- the remaining model errors on Earth: tropics ~3 °C too warm, interior winters too mild, sea-ice extent.

## Addendum (2026-11-22, from E11): snow on high ground, and ice sheets

The ice layer (E11) found that the model decided snow cover from the sea-level temperature, so high ground and
ice-sheet plateaux lost their snow in summer. On Earth that put the Antarctic plateau at −5 °C in January (real
~ −28 °C). `ebm.run` now takes two optional inputs:
- `land_elev`: snow lies where the ground is below `snow_T`, using the ground temperature Ts − 6.5 K/km × land height;
- `ice_sheet`: ice sheets keep the ice albedo all year.

The temperature layer passes the land height always and the ice sheets when it is given `--ice`. On Earth's
stations the change is small (RMSE 5.45 / 3.96 / 3.64 → 5.34 / 3.98 / 3.52 °C), and the coefficients were not
refitted. Details and Antarctic and Greenland values are in E11.

The world's climate code moved to `scripts/world_climate.py`, shared with the precipitation and ice layers.

## Addendum 2 (2026-11-22, from E15): is the North's summer a model error?

The vegetation layer (E15) found the North of Westeros forest-tundra and bog, and the Haunted Forest tundra. The
owner had the climate of the North re-examined before accepting that. Nothing in the model was changed.

- **Stations.** The "subarctic summers 4–8 °C too cold" of the first run belonged to the one-layer model. With the
  two-layer model and the station correction, the July error north of 55° N averages +0.1 °C (leave-one-out):
  - maritime stations are right or too warm: Bergen +0.1, Reykjavik −0.5, Helsinki +0.2, Tromsø +1.2, Murmansk +4.4,
    Nuuk +5.3, Iqaluit +6.8 °C;
  - only the North American interior is too cold: Fairbanks −3.9, Yellowknife −3.3 °C.
- **Station correction.** In July it cools the North by 0.5–1.4 °C on average (Winterfell 2.1 °C, where the
  continentality index is 0.57). Fading its continentality term poleward of 45–70° N does not improve the
  leave-one-out error (3.23 → 3.25–3.32 °C), so the data do not support weakening it.
- **The seas.** July sea-level temperatures:

| Latitude | World, sea west of the North | World, sea east of the North | World, land of the North | Earth model, NE Atlantic | Earth model, Europe |
|---|---|---|---|---|---|
| 56° N | 10.0 | 5.2 | 10.1 | 13.0 | 20.3 |
| 60° N | 7.7 | 3.5 | 8.9 | 9.4 | 18.8 |
| 64° N | 6.0 | 0.6 | 7.0 | 7.7 | 17.6 |

  The world's northern sea keeps ice into July (a third of the sea at 67° N, two-thirds and more north of 69° N).
  The North sits between this cold water on both sides, and its land does not warm above the sea as Europe's does.

**Conclusion.** The cool summers of the North are a consequence of the world's geography, not a model error. The owner
kept the climate and extended canon deviation 5 (CONSTRAINTS decision 11) rather than move the map south.

## Sources

- Station climatology: "List of cities by average temperature", Wikipedia (read 2026-11-21).
- Global mean 14 °C for 1951–1980: NASA Earth Observatory, "World of Change: Global Temperatures".
- North, G. R. (1975), J. Atmos. Sci. 32, 2033 — diffusive energy balance models; North & Coakley (1979) — OLR fit.
- Berger, A. (1978), J. Atmos. Sci. 35, 2362 — insolation.
- Semtner, A. J. (1976), J. Phys. Oceanogr. 6, 379 — zero-layer sea ice.
- Stommel, H. (1948), Trans. AGU 29, 202 — wind-driven gyres.
- Q-flux (prescribed ocean heat transport) in slab-ocean models: standard practice in climate modelling (e.g. Hansen
  et al. 1984).
- Sea-ice extent, 1981–2010 averages — NSIDC (Arctic winter 15.5–15.6 M km², ~6.5 M km² at summer's end; Antarctic
  18.5–18.6 and ~2.5 M km²); NASA Earth Observatory.
