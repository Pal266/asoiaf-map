# E01 — Supercontinent before the breakup

Script `scripts/supercontinent.py`, config `data/tectonics.json` (block `supercontinent`), grid 6 × 200 × 200 (~50 km).
Run: `python3 scripts/supercontinent.py OUT_DIR` (~4 s). Output `supercontinent.npz` and a preview image, not stored
(RULES.md rule 4).

## Method

1. **Landmass.** Seven caps (radius 22–38°) scattered within 45° of (0°, 0°), plus large- and small-scale noise.
   The threshold is found by bisection so that the area is 41 % of the surface (Cogley 1984). Only the largest piece
   is kept, and enclosed seas count as continental crust.
2. **Terranes.** 45 seeds placed at random by area. Each cell joins the nearest seed after its position is warped by
   a smooth noise field (±5°), so the boundaries are not straight. Broken pieces join the neighbour they touch most.
3. **Blocks.** The terranes are clustered on the sphere into 11 blocks, the continents that collided to make the
   supercontinent (Pangaea had about ten). The central terranes of a block, 35 % of its area, form the craton,
   dated 1800–3200 Ma. The other terranes are 350–1200 Ma.
4. **Sutures.**
   - Inside a block: the accretion of a terrane is dated at 55–90 % of the age of the younger terrane, and no
     younger than 450 Ma.
   - Between blocks: the blocks join along a spanning tree of their contacts, longest contacts first. The welds are
     dated in random order from 420 Ma down to the last weld at 220 Ma, 120 Myr before the breakup, as Pangaea
     (320 → 200 Ma). Every contact between two growing clusters closes at the collision that joins them.
   - A suture belt is 250 km wide.

## Result (seed 20261122)

| Check | Value | Earth |
|---|---|---|
| Continental crust | 209 M km², 0.410 of the surface | 210.4 M km², 0.41 (Cogley 1984) |
| Blocks | 11, 6–35 M km² | Pangaea: ~10 major blocks |
| Assembly welds | 410, 406, 371, 361, 331, 321, 285, 283, 252, 220 Ma | Pangaea: Caledonian ~420 → Variscan/Alleghanian ~320 → Cimmerian ~220 Ma |
| Cratonic share | 0.45 of continental crust | Precambrian basement underlies roughly half the continents |
| Suture belts | 0.33 of continental crust; assembly sutures 0.14 | — |

The landmass spans from pole to pole, as Pangaea did (~85° N to 90° S), and centres on the equator.

## Open points

- The shape is a compact blob. Pangaea was C-shaped around the Tethys embayment. If the tectonics needs it, a
  tapered embayment can be added in a later run.
- The ages of the cratons and terranes are drawn at random. Later layers read only the age of the last event
  (`event_age_ma`), so it is enough for them.
- The cratonic share is an order-of-magnitude check, not a calibration.

Next: E02, the breakup along the assembly sutures, and plate kinematics.
