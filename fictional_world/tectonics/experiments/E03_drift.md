# E03 — 100 Myr of drift

Script `scripts/drift.py` (reads E01 and E02), config block `drift` in `data/tectonics.json`.
Run: `python3 scripts/drift.py WORK_DIR` (~70 s). Output `drift.npz`, snapshots every 10 Myr and two preview images,
not stored.

## Method

Rigid-plate advection, semi-Lagrangian, on the 6 × 200 × 200 grid in 1 Myr steps from 100 Ma to 0. Each step, every
cell looks back along each plate's rotation; a plate claims the cell if the source cell belonged to it.
- **No claim:** a gap between diverging plates becomes new ocean floor, dated to that step. It joins the plate of
  its neighbours, so a ridge is shared between the two plates.
- **Several claims:** continental crust overrides ocean floor, and younger floor overrides older (denser) floor.
  Between continents, the larger plate overrides.

Events come from the kinematics of the boundaries each step (normal relative speed below −0.5 cm/yr). They do not
come from the claims, because nearest-cell sampling also duplicates cells along transform faults.
- **Continent over ocean:** a continental arc within 350 km of the boundary (Syracuse & Abers 2006: arcs stand
  ~100–300 km behind the trench).
- **Ocean over ocean:** an island arc 100–220 km behind the boundary, on the overriding ocean floor, but only
  where the arc stays at least 5 Myr. The first version marked every cell a migrating trench swept over, and the
  swaths grew up to ~500 km wide (10.4 % of the ocean floor); now they cover 4.0 %.
- **Continent against continent:** a collision orogen within 300 km. Once ~0.5 M km² of crust has been consumed
  (edge length × convergence × time), the two plates weld. A collision removes the slab pull between them, so the
  welded plate takes a new Euler pole by the reorganisation rule, as the India–Asia collision reorganised plate
  motions. Overridden continental crust is not deleted; it thickens the crust in the band (E04).
- **Accretion.** Ocean crust thicker than 20 km (island arcs, oceanic plateaus) is too buoyant to subduct. Cloos
  (1993) puts the limit at ~17 km for plateaus and arcs. When it meets a continent it is accreted, as Japan, the
  Cordilleran terranes and Wrangellia were. It becomes continental crust in front of the continent, at the ocean
  cells nearest to where it arrived (within 400 km), with its thickness and an arc event. Over 100 Myr: 1.4 M km².
- **Post-collisional convergence.** After the weld the collision zone keeps converging inside the welded plate,
  which rigid plates cannot show. India–Asia has gone on since ~50 Ma and is still active; the Alps ran ~35 Myr
  and the Zagros ~25 Myr.
  - For 40 Myr after a weld, the zone stays an active orogen (event updated each step).
  - Its crust thickens by 0.4 km/Myr, up to 72 km: Tibetan crust doubled to ~70 km.
  - The zone widens outwards by 10 km/Myr, as the Tibetan plateau grew at ~1 cm/yr.
  - These values are order-of-magnitude. Shortening conserves volume: the welded plate loses the equivalent area
    at its coast, 3.4 M km² over 100 Myr.
- **Rifted margin:** continental crust within 150 km of new ocean floor, first time only.
- Every continental cell keeps its terrane and the age and kind of its last event (codes in the config). Every
  oceanic cell keeps its birth age.

**Plate motions.** One rule, `breakup.euler_poles`, is used at the breakup and at the reorganisations at 70 and
40 Ma. Earth had global reorganisations at ~100, ~83 and ~50 Ma (Matthews et al. 2012).
- Plates with continental crust move away from the continental crust of all other plates (inverse-square
  repulsion). After the breakup this points away from the supercontinent; later it spreads the continents over the
  sphere.
- Oceanic plates move towards continental crust, where they are subducted.
- Speeds follow E02.

**Area conservation.** Nearest-cell resampling does not conserve area: a rotated grid skips some source cells and
repeats others, and at plate edges the errors do not cancel. In the first version continental crust shrank by 7 %:
8.7 M km² was overridden in collisions and 6.4 M km² was lost to resampling. Now, each step, every plate's
continental area is restored to its value before the step minus the crust overridden in collisions, a conservative
correction as in mass-conserving semi-Lagrangian schemes:
- missing area goes to the plate's own coast, filling the cells with the most continental neighbours first;
- surplus area comes from the most exposed coastal cells, sparing arcs active in the last 5 Myr. Without this,
  a correction made the active margin retreat, and arc belts of different ages lay side by side in stripes.

Over 100 Myr the correction added 17.8 and removed 9.6 M km².

The first version of E03 used poles that stayed constant for 100 Myr and pointed radially away from the centre. It
pushed the continents towards both poles and into each other. With the reorganisations the continents spread
instead.

## Result (seed 20261122), present day

| Check | Model | Earth |
|---|---|---|
| Mean ocean floor age | 50 Myr | 64 Myr today; younger in the Cretaceous after Pangaea's breakup (Müller et al. 2008) |
| Floor 0–20 / 20–40 / 40–60 / 60–80 / 80–100 / >100 Myr | 22 / 23 / 19 / 17 / 14 / 5 % | the area falls with age; ~15 % older than 100 Myr today |
| Continental crust | 0.397 of the surface: −1.5 % overridden in collisions, −1.6 % post-collisional shortening, +0.7 % accreted; volume 0.98 | Cenozoic collisions consumed ~2 % (India–Asia ~2.5 M km², Arabia, the Alps) |
| Plates | 10 (6 welds: 89, 83, 72, 32, 10, 7 Ma); 2 microcontinents survive (0.4 and 3.9 M km²) | — |

The continents spread into several landmasses (E04). The leading edges of the dispersing continents carry
continental arcs, like the Andes and the Cordillera. Island arcs ring the old ocean, and the thickest of them
accrete where they reach a continent.

## Open points

- Young orogens: without post-collisional convergence, the zones froze at the weld and young orogens covered only
  1.4 % of the continental crust. With it they cover 8 %, plus 1 % plateau (E04).
- Of three microcontinents, one welded back after 90 Myr (at 10 Ma). Real microcontinents also sometimes re-accrete.
- Continental crust includes the shelves. How many separate continents are above sea level is decided in E04, with
  sea level and thickness.
- The absolute reference frame (where the spin axis lies) is free: rigid plate motions do not fix it, and true polar
  wander can shift it. E03 keeps the frame of E01; E04 moves the axis (decision 7).

Next: E04, the output on the 10 km grid: crustal thickness, relief of the ocean floor (GDH1), sea level, land, and
provinces in the codes of the later layers.
