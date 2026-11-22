# E04 — The tectonics layer on the 10 km grid

Script `scripts/output.py` (reads E03), config block `output` in `data/tectonics.json`. Run:
`python3 scripts/output.py WORK_DIR` (~25 s).

E04 also changed E03: overridden continental crust is no longer deleted. Its volume thickens the crust in the
collision band, up to 72 km, and the excess counts as eroded.
- Arcs grow by 0.12 km/Myr: ~30 km³ per km of arc per Myr (Reymer & Schubert 1984) over a ~250 km wide arc. Caps:
  65 km on continents, 30 km on island arcs.
- Rifted margins taper to 12 km at the continent–ocean boundary.
- Starting thickness: cratons 40 km, other basement 37, old sutures 40, ocean 7 (Christensen & Mooney 1995).

## Method

1. **Provinces** (coarse grid), in the codes of the first world's history layer:
   - basement of cratons: shield (40 % of cratons, by a noise field) or platform;
   - other basement younger than 1000 Ma and the sutures of the supercontinent: old orogen;
   - rifted margins: 6;
   - arcs and collisions: 5 and 4 when younger than 66 Ma, otherwise old orogen; a young orogen with crust ≥ 60 km
     thick becomes a plateau (7);
   - old basement 100–400 km from young orogens and arcs: foreland basin (8);
   - island arcs on the ocean floor: 5.
2. **Fine grid.** Each 10 km cell takes the coarse cell under its position warped by a smooth fractal field
   (±40 km, wavelengths 120–1000 km). Thickness and ocean floor age are then smoothed within each crust type
   (6 box passes, σ ≈ 20 km).
3. **Ocean floor depth** from age: GDH1 (Stein & Stein 1992).
4. **First isostatic surface:**
   - Continental relief is Airy on crustal thickness (ρc 2830, ρm 3300 kg/m³), shifted so that the mean height of
     continental crust, shelves included, is 120 m (Cogley 1984).
   - Island arcs stand above the floor by the same rule.
5. **Spin axis (decision 7).** Rigid plate motions do not fix the spin axis, and true polar wander can move it
   (Earth: episodes of up to ~1°/Myr, typically 10–20° over 100 Myr). The whole planet is rotated so that the least
   land, by the coarse isostatic surface, lies poleward of 60°, with the axis moving at most 15°. Of 40,000
   candidate axes, those within 0.5 % of the minimum are kept, and the one nearest the simulation's axis is chosen.
   The longitude origin goes to the centroid of the land.

   | Largest polar wander allowed | 0° | 10° | 15° | 30° | 35° | 45° |
   |---|---|---|---|---|---|---|
   | Land poleward of 60° | 11.4 % | 10.3 % | 10.0 % | 9.4 % | 7.0 % | 4.1 % |

   Polar land falls slowly, because large landmasses reach the poles. The first version had no limit, needed 43°
   and reached 4.3 %; the owner chose a smaller turn. Now the axis moves 8.0° (to 82.0° N, 89.2° W of the
   simulation frame), and polar land falls from 14.2 % to 11.1 % of the land (Earth ~17 %).
   - The rotation is the smallest that moves the axis, so the map stays close to the simulation frame. Centring
     the longitudes on the land's centroid made the map jump between runs. The rotation matrix and the Euler poles in
   the new frame are saved in `tectonics.npz`. E01–E03 previews stay in the simulation frame.
6. **Sea level** is 0.69 × (Earth's mean basin depth − this world's), the water-loading factor of Müller et al.
   (2008). Earth's ocean floor area falls about linearly with age to ~180 Myr (Parsons 1982).

## Output (`tectonics.npz`, 6 × 1000 × 1000, not stored)

- `crust` (0 ocean, 1 continental);
- `province` (codes 0–9);
- `event_age_ma` (last event on continents and island arcs, floor age elsewhere);
- `ocean_depth_m` (below this world's sea level; NaN on continental crust);
- `thickness_km`;
- `plate`;
- `isostatic_elevation_m`;
- `sea_level_m`;
- `frame_rotation`;
- `euler_rad_per_myr`;
- from E05: `hotspot_track_ma`, `lip`, `hotspot_swell_m`, and `active_hotspots.json`.

`volcanoes.npz` lists the hotspot volcano centres (E05) and the volcanic fronts of today's subduction zones: points
150 km inside the overriding plate from each convergent edge (Syracuse & Abers 2006), one per 70 km (arc volcano
spacing ~50–70 km, de Bremond d'Ars et al. 1995). 779 arc volcanoes, 506 of them on continents.

`boundaries.npz` holds the present plate boundaries: position, plates, kind, normal and tangential speed, and the
crust on each side, for the later geology layer.

The whole chain E01 → E04 runs in ~2 minutes, and two runs give identical output.

## Result (seed 20261122)

| Check | Model | Earth |
|---|---|---|
| Continental crust | 0.398 of the surface | 0.41 |
| Mean basement depth of the oceans | 4754 m | 4843 m by the same model with Earth's age distribution |
| Sea level | +60 m above the modern datum | mid-Cretaceous +170 m (Müller et al. 2008); this world's ocean is only ~10 Myr younger than Earth's |
| Land | 0.242 of the surface | 0.29 today |
| Flooded continental crust (shelves, inland seas) | 39 % | 30.6 % today (Cogley 1984) |
| Mean land / ocean depth | 605 m / 3910 m | 840 m / 3680 m |
| Landmasses ≥ 0.5 M km² | 37.1, 23.9, 22.2, 19.5, 11.0, 2.6, 2.1, 1.0, 0.9, 0.6 M km² | 85 (Afro-Eurasia), 42 (Americas), 14 (Antarctica), 8 (Australia) … |

Provinces (share of continental crust): shield 6 %, platform 23 %, old orogen 33 %, young orogen 8 %, arc 15 %,
rifted margin 7 %, plateau 1 %, foreland basin 7 %. Island arcs: 3.1 % of the ocean floor. Plumes, oceanic plateaus and flood basalts: E05.

The numbers are from the final run, with E05 and every E03 fix.

## Open points

- Land stands low (mean 604 m) because the surface has only crustal isostasy. The relief layer adds mountain
  relief, erosion and the thermal uplift of young margins.
- Wide inland seas flood the central continent, much like the Cretaceous seas over Africa and North America.
- Polar land: 11.1 % of the land poleward of 60° after decision 7, less than Earth's ~17 %.
- Island-arc swaths: fixed in E03 (3 % of the ocean floor instead of 10 %).
- Plumes, oceanic plateaus and hotspot tracks: added in E05.
