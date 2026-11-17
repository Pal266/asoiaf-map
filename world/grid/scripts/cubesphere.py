"""Equiangular cubed-sphere grid (world/grid/GRID.md).

Layout
------
Six faces, each N x N cells; arrays are indexed [face, j, i].
On a face, the equiangular coordinates are a (along i) and b (along j), both in [-pi/4, pi/4]; t = tan(a), u = tan(b).
The point on the cube, before normalisation to the sphere (x towards lon 0 / lat 0, y towards lon 90 E, z north):

  face 0  centre lon   0, lat   0   ( 1,  t,  u)   i east, j north
  face 1  centre lon  90, lat   0   (-t,  1,  u)   i east, j north
  face 2  centre lon 180, lat   0   (-1, -t,  u)   i east, j north
  face 3  centre lon -90, lat   0   ( t, -1,  u)   i east, j north
  face 4  north pole                (-u,  t,  1)   i towards lon 90 E, j towards lon 180
  face 5  south pole                ( u,  t, -1)   i towards lon 90 E, j towards lon 0

Cell (i, j) spans a in [-pi/4 + i d, -pi/4 + (i + 1) d], d = (pi/2) / N, and the same for b with j.
Cell areas are exact: the solid angle of the gnomonic rectangle from (0, 0) to (t, u) is
atan(t u / sqrt(1 + t^2 + u^2)), and a cell is the difference of four such corners.
"""

import numpy as np

R_KM = 6371.0          # world/planet/PLANET.md
N = 1000               # world/grid/GRID.md


def face_xyz(face, t, u):
    one = np.ones_like(t)
    if face == 0:
        v = (one, t, u)
    elif face == 1:
        v = (-t, one, u)
    elif face == 2:
        v = (-one, -t, u)
    elif face == 3:
        v = (t, -one, u)
    elif face == 4:
        v = (-u, t, one)
    elif face == 5:
        v = (u, t, -one)
    else:
        raise ValueError(face)
    x, y, z = v
    r = np.sqrt(x * x + y * y + z * z)
    return x / r, y / r, z / r


def xyz_to_latlon(x, y, z):
    """Radians."""
    return np.arcsin(np.clip(z, -1, 1)), np.arctan2(y, x)


def edges(n=N):
    """Equiangular cell edges a_0..a_n (radians)."""
    return np.linspace(-np.pi / 4, np.pi / 4, n + 1)


def centres(n=N):
    e = edges(n)
    return 0.5 * (e[:-1] + e[1:])


def cell_latlon(face, n=N, sub=None):
    """Lat/lon (radians) of cell centres of one face, shape (n, n) as [j, i].
    With sub = (si, sj) in [0, 1): a point inside each cell at that fraction of the cell instead of the centre."""
    e = edges(n)
    d = e[1] - e[0]
    fi, fj = (0.5, 0.5) if sub is None else sub
    a = e[:-1] + fi * d
    b = e[:-1] + fj * d
    t, u = np.meshgrid(np.tan(a), np.tan(b))
    return xyz_to_latlon(*face_xyz(face, t, u))


def latlon_to_cell(lat, lon, n=N):
    """Face and cell indices (face, j, i) of points given in radians."""
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    ax, ay, az = np.abs(x), np.abs(y), np.abs(z)
    face = np.where((ax >= ay) & (ax >= az), np.where(x > 0, 0, 2),
                    np.where(ay >= az, np.where(y > 0, 1, 3), np.where(z > 0, 4, 5)))
    t = np.empty_like(x)
    u = np.empty_like(x)
    with np.errstate(divide="ignore", invalid="ignore"):
        pairs = {0: (y / x, z / x), 1: (-x / y, z / y), 2: (y / x, -z / x),
                 3: (-x / y, -z / y), 4: (y / z, -x / z), 5: (-y / z, -x / z)}
    for f, (tt, uu) in pairs.items():
        m = face == f
        t[m] = tt[m]
        u[m] = uu[m]
    d = (np.pi / 2) / n
    i = np.clip(((np.arctan(t) + np.pi / 4) / d).astype(int), 0, n - 1)
    j = np.clip(((np.arctan(u) + np.pi / 4) / d).astype(int), 0, n - 1)
    return face, j, i


def cell_area_km2(n=N, r=R_KM):
    """Area of each cell, shape (n, n) as [j, i]; the same on every face."""
    t = np.tan(edges(n))
    T, U = np.meshgrid(t, t)
    w = np.arctan(T * U / np.sqrt(1 + T * T + U * U))
    omega = w[1:, 1:] - w[1:, :-1] - w[:-1, 1:] + w[:-1, :-1]
    return omega * r * r


if __name__ == "__main__":
    area = cell_area_km2()
    total = 6 * area.sum()
    sphere = 4 * np.pi * R_KM ** 2
    side = np.sqrt(area)
    print(f"cells: 6 x {N} x {N} = {6 * N * N:,}")
    print(f"total area {total:,.0f} km2, sphere {sphere:,.0f} km2, rel. error {abs(total - sphere) / sphere:.2e}")
    print(f"cell side: centre {side[N // 2, N // 2]:.2f} km, corner {side[0, 0]:.2f} km, "
          f"area ratio max/min {area.max() / area.min():.3f}")
