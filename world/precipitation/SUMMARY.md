# Precipitation — summary

- **What it is:** monthly precipitation of every 10 km cell, from a physical model driven by the temperature layer
  and the relief.
- **The model:** a monthly water-vapour budget on the sphere. Water evaporates from the sea and land, is carried by
  the winds and by mid-latitude eddies, and falls where air rises: at the ITCZ, on windward slopes, in monsoon
  inflow. It falls less where air sinks: under the subtropical highs. A Gill circulation supplies the monsoons.
  Water is conserved: global rain equals global evaporation.
- **Earth test:** 109 stations; typical error ×2.5, 56 % within a factor 2. Zonal means within ~20 % of GPCP.
  Known weaknesses: monsoon rainfall 3–8× too low; arid places ~1.7× too wet; some windward west coasts too wet in
  winter.
- **The world:**

  | Place | Year | January | July |
  |---|---|---|---|
  | Winterfell | 690 mm | 66 | 39 |
  | King's Landing | 286 mm | 26 | 16 |
  | Oldtown | 340 mm | 85 | 3 |
  | Sunspear | 221 mm | 36 | 7 |
  | Asshai | 1317 mm | 35 | 126 |

  The canon deserts are arid or semi-arid: Red Waste ~85 mm, Grey Waste ~130, Great Sand Sea ~170, Dornish desert
  ~270 (likely below 200 after the model's Earth bias). The south of Westeros and the Free Cities have wet winters
  and dry summers.
- **Details:** `experiments/E10_precipitation.md`.
