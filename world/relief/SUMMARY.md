# Relief — summary

- **What it is:** bedrock elevation (cell mean) and expected relief inside each cell, for the whole planet on the
  10 km grid. It is built from the plates, canon relief and geological history layers.
- **How it is built:**
  - Plains take their level from the tectonic province.
  - Every range has an Earth analogue: Frostfangs = Scandinavian Mountains, Red Mountains = Pyrenees,
    Bone Mountains = Tian Shan, the southern continent's west coast = Andes, and so on.
  - The sea floor follows the age of the crust (GDH1), with sediment, trenches, shelves, slopes, oceanic
    plateaus and dynamic topography.
  - All analogue summit heights are checked against references.
- **Result:**
  - Mean land elevation 723 m (Earth 840 m with its ice sheets).
  - Shelf area 26.9 M km² (Earth 28.4 M).
  - Mean ocean depth 3.88 km (Earth 3.69 km). Shelf, slopes and the deep tail match Earth's hypsometry within
    1–3 points. The sea-floor sediment matches GlobSed (margins 3044 m, open ocean 404 m, 337 M km³) and is stored
    as its own layer.
  - The canon coastline is kept exactly.
- **Details:** `experiments/E08_relief.md`.
