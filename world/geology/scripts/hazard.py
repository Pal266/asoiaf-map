"""Relative seismic hazard from active plate boundaries and deforming zones (experiment E16).

Earthquakes concentrate on plate boundaries and in the broad zones of active deformation around them (orogens,
rifts, arcs); shaking falls off with distance from the source. The hazard of a place is the strongest of:
- each class of plate boundary, weight w_class x f(v) x exp(-d / L), with d the distance to the nearest boundary of
  that class moving at relative velocity v; f(v) = min(1, v / v_ref) ** beta, since the seismic moment released on a
  fault grows with its slip rate;
- the deforming zones: weight w_zone inside, falling off as exp(-d / L) outside.
Boundary classes are those of Bird (2003): SUB subduction, OCB oceanic and CCB continental convergent, CTF
continental and OTF oceanic transform, CRB continental rift, OSR oceanic spreading ridge. The index (0-1) is cut into
the four classes of the GSHAP map: low, moderate, high, very high.
"""

import numpy as np
from scipy.spatial import cKDTree

R_KM = 6371.0
CLASSES = ("SUB", "OCB", "CCB", "CTF", "OTF", "CRB", "OSR")


def unit(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def distance_km(src_lat, src_lon, lat, lon, max_km=np.inf):
    """Great-circle distance from every query point to the nearest source point, km; inf when there is no source
    within max_km (searching only that far is much faster)."""
    if len(src_lat) == 0:
        return np.full(np.shape(lat), np.inf)
    bound = 2 * np.sin(min(max_km / R_KM, np.pi) / 2) if np.isfinite(max_km) else np.inf
    d, _ = cKDTree(unit(np.asarray(src_lat), np.asarray(src_lon))).query(
        unit(np.ravel(lat), np.ravel(lon)), distance_upper_bound=bound, workers=-1)
    return np.where(np.isfinite(d), 2 * np.arcsin(np.clip(d / 2, 0, 1)) * R_KM, np.inf).reshape(np.shape(lat))


def rate_factor(v, cfg):
    return np.minimum(1.0, np.asarray(v, float) / cfg["rate_ref_mm_yr"]) ** cfg["rate_exponent"]


def index(sources, zone_inside, zone_dist, cfg):
    """sources: list of (class, velocity mm/yr, distance km array); zone_inside: weight of the deforming zone a point
    lies in (0 outside); zone_dist: distance to the nearest deforming zone, km. Returns the hazard index (0-1)."""
    L = cfg["decay_km"]
    h = np.asarray(zone_inside, np.float64).copy()
    for c, v, d in sources:
        h = np.maximum(h, cfg["boundary_weight"][c] * rate_factor(v, cfg) * np.exp(-d / L))
    h = np.maximum(h, cfg["zone_weight_outside"] * np.exp(-zone_dist / L))
    return h


def classes(h, cfg):
    """GSHAP classes 0 low, 1 moderate, 2 high, 3 very high."""
    return np.digitize(h, cfg["class_thresholds"])
