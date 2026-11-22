# Global grid

## Resolution (project owner's decision)

- The layers of the global level are delivered on a single grid with a cell size of **about 10 km** (about 100 km² per
  cell). On a planet with Earth's surface area (about 510 million km²) this is about 5.1 million cells.
- A layer may be **computed** on a coarser grid when the process it models works at a larger scale, and then brought
  to the 10 km grid:
  - tectonic history: about 50 km (plates and mountain belts are hundreds of kilometres across);
  - climate: about 25–50 km (the resolution of global climate models);
  - relief, rivers, erosion and soils: 10 km (coasts, rivers and lakes need it).
- The document of every layer states the grid it was computed on.

## Grid type (project owner's decision)

- **Equiangular cubed sphere.** The sphere is projected onto the 6 faces of a cube, and each face is divided into
  **1000 × 1000** cells at equal angles: 6 000 000 cells. A cell is 10.0 km across at the centre of a face
  (6 371 km × 90° / 1000) and 8.8 km at the corners; the largest cell has 1.41 times the area of the smallest, so
  every computation uses each cell's true area.
- Why: each face is a plain 2D array, so flow, erosion and smoothing algorithms stay simple and fast, and a region is
  easy to cut out for a detailed map. Modern global climate models use this grid (NOAA's FV3 dynamical core).
- Cells on the edges of a face take their neighbours from the adjacent face.
- The numbering and orientation of the faces are fixed in `scripts/cubesphere.py` (its docstring): faces 0–3 are
  centred on the equator at longitudes 0°, 90° E, 180° and 90° W, face 4 is the north pole, face 5 the south pole;
  arrays are indexed `[face, j, i]`. Cell areas are computed exactly (solid angle of each cell); the 6 faces add up
  to the sphere's area to a relative error of 5 × 10⁻¹⁶ (`python3 scripts/cubesphere.py`).
- Neighbours across face seams: `halo_index(n, p)` gives, for each face padded by `p` cells, the real cell under
  every padding cell (equiangular coordinates continued past the edge). `box_smooth` (repeated 3 × 3 average,
  Gaussian with σ ≈ √(2k/3) cells after k passes) and `boundary_cells` (cells of a mask with a neighbour outside it)
  use it, so smoothing and edge detection have no seams.
