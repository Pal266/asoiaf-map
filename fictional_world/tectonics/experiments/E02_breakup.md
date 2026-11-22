# E02 — Plates at the breakup (100 Ma)

Script `scripts/breakup.py` (reads the E01 output), config block `breakup` in `data/tectonics.json`.
Run: `python3 scripts/breakup.py WORK_DIR` (~5 s). Output `breakup.npz` and a preview image, not stored.

## Method

1. **Continental plates.** The supercontinent splits along the assembly sutures between its 11 blocks. Shorter
   sutures reactivate first (Erdős et al. 2025), so the longest block contacts stay welded until 9 continental
   plates remain. Internal sutures of the blocks stay closed.
   - **Microcontinents.** Small non-cratonic terranes on a rifted edge break away as their own plates, as
     Madagascar, the Seychelles and the Lord Howe Rise left Pangaea's fragments. Candidates are under 4 M km²,
     one per plate, and never the last terrane of a plate.
   - Of up to 9 candidates, the 3 that open fastest from the plate they left are kept, after a first pass of
     poles; a sliver pushed back would weld again at once. The poles are then recomputed.
   - Kept: 1.45, 0.54 and 3.79 M km², opening at 2.0, 1.0 and 0.9 cm/yr.
2. **Oceanic plates.** The exterior ocean is shared among 7 plates. Their seeds are at least 20° apart, and the
   boundaries are warped Voronoi edges.
3. **Euler poles.**
   - Each continental plate moves away from the supercontinent's centre, at 1.5–6 cm/yr. Plates further out move
     faster; the speed drops by up to 40 % with the cratonic share, after Zahirovic et al. (2015): ~2.8 cm/yr for
     plates more than 25 % cratonic, ~5.8 for plates with any craton.
   - Each oceanic plate moves towards the supercontinent at 5–9 cm/yr, so the exterior ocean is consumed at its
     margins. Such plates are bounded by subduction zones, and Zahirovic et al. give them ~8.5 cm/yr.
   - Every plate gets a random spin about its centroid (σ 0.12°/Myr).
4. **Boundaries.** The relative velocity across each cell edge is split into a normal and a tangential part.
   - Divergent: normal speed > 0.5 cm/yr.
   - Convergent: normal speed < −0.5 cm/yr.
   - Transform: anything else, or a normal speed less than half the tangential.
5. **Exterior ocean floor.** The age of a cell is its distance to the nearest divergent ocean–ocean boundary divided
   by that boundary's half spreading rate, capped at 200 Myr. Older floor would already have been subducted, and
   the Panthalassa floor of the Jurassic reached ~200 Myr.

## Result (seed 20261122)

| Check | Model | Earth |
|---|---|---|
| Plates | 12 continental (3 of them microcontinents) + 7 oceanic | 16–46 over the last 410 Myr; fewer at a supercontinent stage (Matthews et al. 2016) |
| Continental plates, median RMS | 4.6 cm/yr | cratonic plates ~2.8, any craton ~5.8 (Zahirovic et al. 2015) |
| Oceanic plates, median RMS | 6.5 cm/yr | fastest ~8.5 |
| All plates, area-weighted RMS | 5.6 cm/yr | Mesozoic–Cenozoic ~6, median ~4 (Young et al. 2019; Zahirovic et al. 2015) |
| Rifts (continent–continent) | 62 % divergent, 31 % transform, 7 % convergent | young rifts open obliquely, with transform segments (Atlantic) |
| Continental margins | 49 % convergent, 34 % transform, 16 % divergent | a supercontinent is ringed by subduction |
| Exterior ocean floor | mean 80 Myr, median 56; 14 % at the 200 Myr cap | — |

The central plate (8 % of the surface, at 10° from the centre) barely moves, like Africa after Pangaea. The others
spread outwards. The oldest floor lies along the continental margins, where it is subducted. The ridges of the
exterior ocean lie on the far side of the planet, as the Pacific ridges did around Pangaea.

## Open points

- The numbers above are from the final run (with the repulsion rule of E03 and microcontinents). 7 % of the rift
  edges converge: a short section of rift closes. E03 will show whether this matters; if so, it can
  become a small fold belt.
- 13 % of the continental margins diverge: there the oceanic plate moves away from the continent, which would open
  a new ridge near the coast.
- The poles are constant for all 100 Myr. Real plates change motion every ~10–50 Myr; a reorganisation can be added
  in E03 if the drift looks too regular.

Next: E03, 100 Myr of drift in 1 Myr steps.
