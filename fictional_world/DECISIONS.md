# Owner's decisions

This world has no canon. The owner's design decisions play its part: every later layer must respect them, and a
layer that cannot (Rule 1: science first) stops for the owner's choice (Rule 3).

1. **Folder** (2026-11-22): the world lives in `fictional_world/`, on branch `fictional_world`.
2. **A planet like Earth** (2026-11-22): Earth's size, orbit, axial tilt, atmosphere, age and Moon
   (`planet/PLANET.md`). The climate and other models can then be tested on Earth without recalibration.
3. **Recent breakup of a supercontinent** (2026-11-22): the world is seen some tens of millions of years after a
   supercontinent began to break up: many continents, long coasts, inland seas and rifts.
4. **The game region is chosen after the climate** (2026-11-22): the world is built through the climate and the
   hydrology first, and the region is picked among real candidates.
5. **The global grid of the first world** (2026-11-22, taken over by default; the owner may change it):
   equiangular cubed sphere, 6 × 1000 × 1000 cells, ~10 km (`grid/GRID.md`).
6. **The breakup began 100 Myr ago** (2026-11-22, research `tectonics/research/R01_breakup_age.md`): the world
   matches Pangaea in the mid-Cretaceous, ~80–100 Myr after its breakup began. The tectonic simulation runs from an
   assembled supercontinent through 100 Myr of rifting and drift.
7. **The spin axis is moved a little towards less polar land** (2026-11-22, owner's choice after tectonics E04):
   rigid plate motions do not fix the spin axis, and true polar wander can move it. The axis may move at most 15°,
   a typical true polar wander over 100 Myr. Within that limit it goes where the least land lies poleward of 60°.
   It moved 8.0°, and polar land fell from 14.2 % to 11.1 % of the land (Earth ~17 %). A first version without the
   limit reached 4.3 % but needed 43°; the owner found that too steep.
