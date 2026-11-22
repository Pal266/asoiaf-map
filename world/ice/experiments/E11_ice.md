# E11 — Ice sheets and glaciers in equilibrium with the climate

Status: **done** — run j is the layer (owner decision 2, canon deviation 6).

## Goal

Where the world's climate keeps permanent ice, and how thick it is. Canon checks: the Land of Always Winter is
permanently frozen (WE7); the haunted forest, the North and Winterfell are not under ice.

## Model

**Surface mass balance** (`scripts/smb.py`): the positive-degree-day model of ice-sheet modelling.
- Snowfall: the share of monthly precipitation that falls as snow, all snow below −1 °C, all rain above 3 °C.
- Melt: positive degree days from the monthly mean and a daily standard deviation σ (Calov & Greve 2005), at
  ddf_snow for snow, then ddf_ice for ice.
- Refreezing: up to 60 % of the year's snowfall (Reeh 1991).
- A raised (ice) surface is colder by 6.5 K/km and drier by 7 %/K (Huybrechts 2002).

**Glaciers** outside the sheets: the share of a cell above its equilibrium line, from its sub-grid relief, divided
by an accumulation-area ratio of 0.6.

**Ice sheets** (`scripts/ice.py`):
- connected areas of positive balance of at least 5000 km²;
- grown into the neighbouring ablation zone, least negative first, while the sheet's total balance (surface balance
  minus calving where it reaches the sea) stays positive;
- calving: speed ∝ water depth (Brown, Meier & Post 1982), c = 2.5 /yr, through the margin thickness;
- thickness: perfectly plastic ice (Nye 1951) with yield stress τ0, integrated inwards from the margin over the
  bed, with local isostatic depression of the bed (ρ_ice / ρ_mantle).

**Coupling with the climate.** Ice changes the climate (albedo, height) and the climate changes the ice. The layer
iterates from bare ground:
1. climate of the current surface: the E09 temperature model with the land height and the ice sheets' albedo, and
   the E10 moisture model (`world/climate/scripts/world_climate.py`);
2. inner loop: with that climate, ice and its own height to equilibrium;
3. repeat until the sheets change by less than 100 000 km² between passes.

### Changes to the temperature model (E09)

The EBM decided snow cover from the sea-level temperature, so high ground and ice-sheet plateaux were dark in
summer. On Earth that left the Antarctic plateau at −5 °C in January (real ~ −28 °C) and no ice in Antarctica's
interior. Two optional inputs were added to `ebm.run`:
- `land_elev`: snow lies where the ground, not sea level, is below snow_T;
- `ice_sheet`: ice sheets keep the ice albedo all year.

| Earth test (EBM, 78 stations) | RMSE Jan / Jul / year °C | Global mean | South Pole Jan / Jul | Summit Jan / Jul |
|---|---|---|---|---|
| before | 5.45 / 3.96 / 3.64 | 17.7 | −6 / −44 | −40 / −6 |
| land height | 5.34 / 3.98 / 3.52 | 16.1 | −17 / −48 | −42 / −9 |
| land height + observed ice | 5.41 / 4.04 / 3.60 | 15.8 | −18 / −48 | −42 / −12 |

Real: South Pole −28 / −60, Summit ~ −40 / −11. The stations barely change. The temperature and precipitation
layers now use the land height. World effect without ice: global mean 15.6 → 15.2 °C, haunted forest's warmest
month 7.6 → 7.4 °C, Winterfell −9.8 / 11.1 °C.

The climate models of the world were moved into one module (`world_climate.py`) used by the temperature,
precipitation and ice layers, so that all three see the same climate for the same surface. The temperature and
precipitation layers take an optional `--ice ice.npz`.

## Earth tests (`scripts/earth_ice_validation.py`)

Inputs: Earth climate from the E09 and E10 models, downscaled to 10′; ETOPO1 surface; observed ice from Natural Earth
1:10m glaciated areas (14.65 M km² on land).

### a. Balance on today's surface

Is the balance positive where ice is today? Area-weighted critical success index (hits / (hits + misses + false
alarms)):

| σ K | ddf_snow mm/K/d | CSI | hit rate | false alarms M km² | Antarctica k km² | Greenland k km² | High Mountain Asia k km² |
|---|---|---|---|---|---|---|---|
| observed | | | | | 12 242 | 1 862 | 145 |
| 3.5 | 3.0 | **0.81** | 0.84 | 0.52 | 11 512 | 1 102 | 202 |
| 4.0 | 3.5 | 0.78 | 0.80 | 0.37 | 10 931 | 1 005 | 140 |
| 5.0 | 3.0 | 0.72 | 0.73 | 0.31 | 9 985 | 926 | 137 |
| 5.0 | 4.0 | 0.70 | 0.71 | 0.22 | 9 654 | 833 | 79 |

Before the EBM change: hit rate 0.31 (Antarctic interior melting).

Small maritime ice caps are missed: Iceland, Svalbard, Patagonia, coastal Alaska, ~0.3 M km². Their precipitation is
3–10× too low, because a 2° moisture model cannot resolve their narrow mountains.

### b. Ice thickness

Plastic ice on the observed masks against the ETOPO1 surface. The bed under the ice is not in the inputs, so 0 and
300 m were tried:

| τ0 kPa | Greenland mean / max surface m | Antarctica mean / max m |
|---|---|---|
| observed | 1920 / 3224 | 2337 / 4087 |
| 50 | 979–1279 / 1981–2281 | 1857–2157 / 3260–3560 |
| **75** | 1199–1499 / 2426–2726 | 2274–2574 / 3992–4292 |
| 100 | 1384–1684 / 2801–3101 | 2626–2926 / 4610–4910 |

τ0 = 75 kPa: Antarctica is right; Greenland, with its high mountain rim, is too low.

### c. Coupled, from bare ground (hysteresis)

The full loop on Earth with all ice removed (σ 3.5):
- **Antarctica regrows:** 9.2–13.5 M km² (observed 12.2);
- **Greenland does not** (0): consistent with models in which an ice-free Greenland does not regrow under
  pre-industrial climate (Ridley et al. 2010);
- Canada, Siberia and Scandinavia stay ice-free;
- **false:** an ice sheet of ~2 M km² on the Tibetan plateau. On the plateau the model is 2–5 °C too cold in summer:
  the free-air lapse rate from sea level ignores the warmth of high plateaux.

Coupled, from a glacial start (ice wherever the bare-ground annual mean is below −5 °C): the ice grew to 25 M km² and
kept growing. Earth's real hysteresis allows large northern ice sheets near present insolation (Abe-Ouchi et al.
2013), so this test does not decide between the starting states.

## World runs (`ice.py`)

| Run | Start | σ K / ddf_snow | Sheets M km² | Ice south to | Notes |
|---|---|---|---|---|---|
| a | bare | 3.5 / 3.0, τ0 100, no calving, no inner loop | 3.5, still growing after 8 passes | — | the south polar continent stays bare |
| g | glacial | 3.5 / 3.0, τ0 100 | 21.3 | 48° N | south polar continent covered; max 7.5 km thick |
| b | bare | 3.5 / 3.0, τ0 100, inner loop | 7.0 | ~50° N | the North and Winterfell under 2–2.5 km of ice |
| c | bare | 3.5 / 3.0, τ0 75, calving | 6.4 | ~50° N | as b |
| **e** | bare | **5.0 / 4.0**, τ0 75, calving | **3.0** | **~64–68° N** | Land of Always Winter covered (99.7 %); the Wall's latitude 37 % ice, 65–66° N 53 %, 66–68° N 86 % |

The south polar continent stays ice-free from bare ground in every run. Its summer is +8 °C at the pole, with
~70–150 mm of precipitation.

**Owner decision 1 (2026-11-22): canon decides within the model's uncertainty** (CONSTRAINTS decision 8): σ 5 K,
ddf_snow 4. In the same pass `ocean_e0_mm_day` of E10 was refitted (E10 addendum).

| Run | Start | σ K / ddf_snow | Sheets M km² | Land of Always Winter | Ice south of the Wall | Notes |
|---|---|---|---|---|---|---|
| f | bare | 5.0 / 4.0 | 3.2 | 99.8 % | 11 % of 60–63.75° N (northern mountains, 47 000 km²) | haunted forest belt 50–67 % ice; south polar continent bare; plateau sheet 0.59 M km² at 40° S |
| h | history | 5.0 / 4.0 | 18.3 | 99.9 % | 40 % of 60–63.75° N | south polar continent covered (14.0 M km², up to 6.6 km); the cold it adds spreads the northern sheet; new sheet 0.78 M km² in north-east Essos |
| i | history | 6.0 / 5.0 | 14.6 | 0 % | 0 | no northern sheet |
| j | history | 5.5 / 4.5 | 14.7 | 2 % | 0 | no northern sheet |

"history" start: the south polar continent starts glaciated, everything else bare. That continent has sat at the pole
inside the Antarctic ring ridge through the whole geological history (E06, E07), as Antarctica has since the
Eocene–Oligocene. From bare ground its summer at the pole is +8 °C and it stays bare. Earth's Antarctica in the
same model is equally warm when bare (+10 °C at sea level in January), yet it regrows in test c, so the two
outcomes sit close to a threshold.

Temperature with the ice (`climate.py --ice`):

| Run | Global mean °C | Castle Black Jan / Jul | Winterfell Jan / Jul | Haunted forest warmest month |
|---|---|---|---|---|
| no ice | 15.2 | −18.5 / 7.9 | −9.8 / 11.1 | 7.4 |
| f | 14.9 | −23.0 / −0.1 | −11.4 / 10.1 | −0.2 |
| h | 13.1 | −24.7 / −2.5 | −12.2 / 9.5 | −1.7 |

Earth balance test (a) for the stronger melt: σ 5.5 / 4.5 → CSI 0.64, Earth's ice 9.7 M km² (−34 %);
σ 6 / 5 → CSI 0.58, 8.7 M km² (−41 %).

**The northern ice is bistable.** Between σ 5 and 5.5 K the Land of Always Winter's sheet jumps from reaching the
Wall to not existing. No parameter value gives a sheet confined to the Land of Always Winter. This is the small
ice-cap instability of energy-balance models (North 1984): a cap smaller than a critical size cannot sustain
itself through its own albedo, and a larger one spreads.

**Owner decision 2 (2026-11-22): no northern ice sheet** — σ 5.5 K, ddf_snow 4.5, history start (run j).
CONSTRAINTS decision 8, canon deviation 6.

## Result (run j)

| Ice sheet | Area | Mean / max thickness | Bed | Comment |
|---|---|---|---|---|
| south polar continent | 13.99 M km² | 2.9 / 6.6 km | 550 m | Antarctica: 14.0 M km², 2.1 / 4.8 km; plastic ice without ice streams is too thick |
| high plateau, 40° S 126° W | 0.63 M km² | 0.8 / 3.1 km | 3.8 km | probably too large: on Earth's Tibet the model is 2–5 °C too cold (test c) |
| north polar island, 85° N 61° E | 59 000 km² | 0.8 / 1.3 km | 240 m | |
| Land of Always Winter mountains, 72° N | 28 000 km² | 0.6 / 0.9 km | 1.3 km | an ice cap |
| tropical plateau, 11° S 116° W | 30 000 km² | 0.6 / 1.1 km | 5.3 km | above 5 km |

- Ice volume 42 M km³ (Earth ~30); 106 m of sea-level equivalent. Mountain glaciers outside the sheets 0.12 M km².
- No ice in Westeros south of 70° N. Climate with this ice: global mean 13.4 °C; Castle Black −19.0 / 7.6 °C,
  Winterfell −11.1 / 10.7 °C (Jan / Jul); haunted forest warmest month 7.0 °C; Land of Always Winter 2.7 °C.
- Convergence: 14.52 → 14.66 → 14.72 → 14.74 M km² in four climate passes (last change 19 000 km²).

## Still open

- The plateau cold bias: raising the station correction's height cap fixes Tibet but warms the tropical Andes by
  4–7 °C (tested caps 2.5, 3, 4 km; kept 2.5).
- Small maritime ice caps are missed, because the precipitation is too coarse for narrow mountains.
- Plastic sheets are too thick in their interiors: no ice streams.

## Sources

- Reeh, N. (1991), Polarforschung 59, 113 — degree-day model, refreezing.
- Calov, R. & Greve, R. (2005), J. Glaciol. 51, 173 — positive degree days with a daily temperature spread.
- Huybrechts, P. & de Wolde, J. (1999), J. Climate 12, 2169; Huybrechts, P. (2002), Quat. Sci. Rev. 21, 203.
- Budd, W. F. & Smith, I. N. (1979) — elevation desertification.
- Meier, M. F. & Post, A. S. (1962), IASH 58 — accumulation-area ratio.
- Nye, J. F. (1951), Proc. R. Soc. A 207, 554 — perfectly plastic ice sheet.
- Brown, C. S., Meier, M. F. & Post, A. (1982), USGS Prof. Paper 1258-C — calving speed and water depth.
- Ridley, J. et al. (2010), Clim. Dyn. 35, 1049 — Greenland does not regrow from ice-free.
- Abe-Ouchi, A. et al. (2013), Nature 500, 190 — northern ice-sheet hysteresis.
- North, G. R. (1984), J. Atmos. Sci. 41, 3390 — small ice-cap instability.
- Natural Earth 1:10m glaciated areas (public domain), github.com/nvkelso/natural-earth-vector.
- ETOPO1 (NOAA), 10′ grid via fatiando-data, DOI 10.5281/zenodo.5882203.
