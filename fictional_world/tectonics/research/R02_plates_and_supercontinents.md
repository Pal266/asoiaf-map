# R02 — Numbers for the plate simulation

## Continental crust

- Continental crust (continents and their shelves) covers ~41 % of Earth's surface: 210.4 million km², in 14
  continents; 30.6 % of it is submerged (Cogley 1984, Rev. Geophys. 22, 101).
- Its area has stayed at ~40 % of the surface for the last 3 billion years (Cawood & Hawkesworth 2019, Gondwana Res.).
- So the supercontinent of this world holds ~41 % of the surface, and its pieces should make roughly a dozen
  continents.

## Assembly and breakup

- Pangaea was assembled from continental blocks by ~320 Ma and began to break up at ~200 Ma: it lasted ~120 Myr as
  a whole (Leleu et al. 2015, Tectonics 34, 1009).
- Rifts and passive margins often follow the old suture zones where continents had collided: the North Atlantic
  opened along Caledonian and Variscan sutures (Wilson 1966; Buiter & Torsvik 2014; Schiffer et al. 2020, Gondwana
  Res.). Not every suture is reused, and the break can follow an older structure. Among several sutures, the oldest,
  shortest and best-oriented one tends to reactivate (Erdős et al. 2025, Solid Earth 16, 841).
- A supercontinent is surrounded by an exterior ocean with old floor and subduction at its margins. The young
  interior oceans open between the dispersing pieces (Murphy, Nance & Cawood 2009).

## Plate speeds

- The median RMS plate speed over the last 200 Myr is ~4 cm/yr. Plates with any cratonic part move at a median
  ~5.8 cm/yr, plates with more than 25 % cratonic area at ~2.8 cm/yr. The fastest plates (~8.5 cm/yr) are mostly
  oceanic and bounded by subduction zones. Plates with more than half continental area do not exceed ~10 cm/yr
  (Zahirovic et al. 2015, EPSL 418, 40).
- Present-day continental RMS speed ~3 cm/yr, all plates ~5 cm/yr. The number of plates over the last 410 Myr is
  16–46, with 2.7–6.6 small plates per large one (Matthews et al. 2016, Global Planet. Change 146, 226).

## A method for simulated plate histories

Procedural plate simulations move the crust of each rigid plate by rotation about its Euler pole on a sphere grid:
- where plates diverge, new ocean floor forms at the ridge;
- where they converge, one side is subducted or the continents collide and fold.

PlaTec (Viitanen 2012) and "Procedural tectonic planets" (Cortial, Peytavie, Galin & Guérin 2019, Computer Graphics
Forum 38) work this way. This layer follows that scheme with Earth's numbers above.
