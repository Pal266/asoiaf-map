# E05 — Mantle plumes: hotspots, tracks, oceanic plateaus and flood basalts

Added to `scripts/drift.py` (E03) and `scripts/output.py` (E04), config block `plumes`. Research: `research/R03_plumes.md`.

## Method

- **Plumes** sit in the simulation (mantle) frame: 7 primary and 90 secondary, at least 8° apart, placed at random.
  Each drifts slowly in a random direction: primary plumes up to 1 cm/yr, secondary ones up to 2 cm/yr. Courtillot
  et al. (2003) find the primary hotspots moving slowly within each hemisphere, and up to 5 cm/yr between the
  hemispheres before 50 Ma.
- **Primary plumes** start with a plume head, which makes a large igneous province: a cap of 0.5–2 M km² (log-uniform).
  - Three start within 10 Myr of the breakup, as LIPs cluster at supercontinent breakup; the other four at random
    times.
  - They then stay active to the present.
  - Oceanic crust under a head thickens to 25–35 km and becomes an oceanic plateau (province 9).
  - Continental crust gets +3 km of flood basalt and underplating (~+430 m by Airy, like the Deccan plateau at
    ~0.5–1 km), and the event "large igneous province" (code 8).
  - The outline is irregular (radius ±30 % by a smooth noise field), and the lavas thin over the outer 250 km. A
    first version with a sharp circular rim left a flat block with straight scarps in the relief.
- **Secondary plumes** live 10–60 Myr, starting at random times. About a third are active at any time, so ~40
  hotspots are active today, as on Earth.
- **Tracks.** Each step, crust within 100 km (primary) or 60 km (secondary) of an active plume records the time.
  - Oceanic crust thickens to 18 or 12 km.
  - Continental crust gets the event "hotspot volcanism" (code 7).
- **Volcano centres.** Every 50 km (primary) or 80 km (secondary) of plate motion over a plume, a volcano centre is
  born and then rides on its plate. Hawaiian volcanoes stand ~40–60 km apart (ten Brink 1991). Centres whose crust
  passes to another plate are dropped (subducted or overridden). 1178 remain today, 306 of them from primary
  plumes. The relief layer builds the edifices on them.
- **Swells** of the plumes active today: dynamic uplift of 1000 m (primary) or 500 m (secondary), Gaussian with a
  half-width of 500 km, in the output frame.

## Result (seed 20261122)

| Check | Model | Earth |
|---|---|---|
| Active hotspots today | 40 (7 primary + 33 secondary) | 49, of which 7 primary (Courtillot et al. 2003) |
| LIPs in 100 Myr | 7: 97, 96, 95 Ma (continental, at the breakup), 63 (oceanic), 49 (oceanic), 28 (continental), 15 Ma (oceanic) | ~1 per 20 Myr, clustered at breakups |
| LIP area | 0.6–2.0 M km² | Ontong Java 1.86 M km² |
| Oceanic plateaus | 0.55 % of the ocean floor (some accreted to continents) | — |
| Hotspot tracks | 2.1 % of the ocean floor | — |
| Flood basalts | 3.0 % of continental crust | — |

The three breakup LIPs fall on the continents, like the Central Atlantic and Karoo provinces at Pangaea's breakup.

## Open points

- **No volcanic islands rise above sea level in this layer.** The isostatic surface spreads a track's thickened
  crust over 10 km cells with Airy compensation. A real volcano is a ~100 km load held up by the flexure of the
  plate, and its summit rises 9–10 km above the floor (Hawaii). The relief layer builds the edifices on the volcano
  centres (`volcanoes.npz`).
- Plumes do not trigger rifting, and plate speeds do not respond to plume heads. Zahirovic et al. (2015) see
  short accelerations at LIP times.
