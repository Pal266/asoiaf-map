# Tectonics — method (planned)

Grid: the cubed sphere at 6 × 200 × 200 cells (~50 km, `../grid/GRID.md`: tectonics may be computed at ~50 km),
brought to the 10 km grid at the end.

1. **Supercontinent (experiment E01).** All continental crust (41 % of the surface) in one landmass, built from
   terranes:
   - cratons (Archean and Proterozoic cores);
   - younger accreted terranes;
   - the suture belts between them, with the ages of the collisions that welded them. The last assembly is
     ~120 Myr before the breakup, as for Pangaea.

   An exterior ocean with old floor surrounds it.
2. **Breakup (E02).** At 100 Ma the supercontinent splits along chosen sutures (Wilson cycle) into ~10–14
   continental plates. The exterior ocean is divided into oceanic plates. Each plate rotates about its own Euler
   pole, at Earth's speeds: continental plates ~3–6 cm/yr, oceanic up to ~9 cm/yr.
3. **Drift, 100 Myr in 1 Myr steps (E03).** The plate motions are recomputed at reorganisations (70 and 40 Ma). The crust moves semi-Lagrangian, each plate by its rotation:
   - gaps become new ocean floor, dated by the step;
   - where plates converge, oceanic floor under a continent is subducted (active margin, arc);
   - older oceanic floor under younger is subducted (island arcs);
   - continent against continent collides (orogen, dated), and the plates weld.
4. **Mantle plumes (E05)**, inside the drift: hotspots fixed in the mantle, tracks, plume heads as large igneous
   provinces, and swells of the active hotspots.
5. **Output (E04).** On the 10 km grid:
   - land and sea from crust and sea level;
   - crust type and ocean floor age and depth (GDH1);
   - plates and their boundaries with relative velocities;
   - tectonic provinces and the age of their last event, in the codes the later layers of the first world read.

Each step is tested against Earth's numbers:
- continental area and number of continents;
- ocean floor age distribution (mean age ~64 Myr today; younger after a breakup);
- plate speeds;
- shares of active and passive margins.
