import math
import pickle

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


def box_to_convex_polygon(box):
    """Convert an lsst.sphgeom.Box to a ConvexPolygon.

    Parameters
    ----------
    box : lsst.sphgeom.Box
        A box on the sphere.

    Returns
    -------
    ConvexPolygon
        The equivalent polygon.
    """
    from lsst.sphgeom import ConvexPolygon, LonLat, UnitVector3d

    if box.isEmpty():
        raise ValueError("Cannot convert an empty Box to a ConvexPolygon.")

    # Get the corners of the box
    lon_a, lon_b = box.getLon().getA().asRadians(), box.getLon().getB().asRadians()
    lon_min = min(lon_a, lon_b)
    lon_max = max(lon_a, lon_b)
    lat_a, lat_b = box.getLat().getA().asRadians(), box.getLat().getB().asRadians()
    lat_min = min(lat_a, lat_b)
    lat_max = max(lat_a, lat_b)
    # todo : this may be an improper assumption, considering RA wrap around!!

    # Convert corners to UnitVector3d
    corners = [
        LonLat.fromRadians(lon_min, lat_min),
        LonLat.fromRadians(lon_max, lat_min),
        LonLat.fromRadians(lon_max, lat_max),
        LonLat.fromRadians(lon_min, lat_max),
    ]
    vertices = [UnitVector3d(corner) for corner in corners]
    return ConvexPolygon(vertices)


def unit_vector3d_to_radec(vec):
    """Convert a UnitVector3d to RA/Dec degrees.

    Parameters
    ----------
    vec : lsst.sphgeom.UnitVector3d
        A 3D unit vector.

    Returns
    -------
    tuple of float
        (RA in degrees [0, 360), Dec in degrees)
    """
    from lsst.sphgeom import LonLat

    lonlat = LonLat(vec)
    ra = lonlat.getLon().asDegrees() % 360.0
    dec = lonlat.getLat().asDegrees()
    return ra, dec


def load_pickle_skymap(path):
    """Load a (raw/unconverted) SkyMap object from a pickle file.

    Parameters
    ----------
    path : str or Path
        Path to the pickle file containing the SkyMap object.

    Returns
    -------
    lsst.skymap.SkyMap
        The loaded SkyMap object.
    """
    with open(path, "rb") as f:
        return pickle.load(f)


def radians_to_degrees(radians):
    """Convert radians to degrees."""
    return radians * (180.0 / math.pi)


class IterateTractAndRing:
    """An iterator to traverse through tract ids and ring ids in a skymap.

    This iterator yields tuples of (tract_index, ring_index) for each tract in the skymap.

    Note that the first tract is the south pole, and is considered to be in ring -1, but is
    assumed to be excluded from ``ring_nums``. If ``add_poles`` is True, the iterator will
    include the south pole (tract 0) and the north pole (last tract) in the iteration.

    Parameters
    ----------
    ring_nums : list of int
        A list where each element represents the number of tracts in each ring (eg, [5, 10, 15] will
        be interpreted as 5 tracts in ring 0, 10 in ring 1, and 15 in ring 2).
    add_poles : bool, optional
        If True, the iterator will include the south pole (tract 0) and the north pole (last tract).
        If False, it will only iterate through the rings. Default is True.
    """

    def __init__(self, ring_nums, add_poles=True):
        self.ring_nums = ring_nums
        if add_poles:
            self.total_tracts = sum(ring_nums) + 2
            self.current_tract = 0
            self.current_ring = -1
        else:
            self.total_tracts = sum(ring_nums)
            self.current_tract = 1
            self.current_ring = 0

    def __iter__(self):
        return self

    def __next__(self):
        # End iteration if we have processed all tracts.
        if self.current_tract >= self.total_tracts:
            raise StopIteration
        tract_and_ring = (self.current_tract, self.current_ring)

        # Increase tract.
        self.current_tract += 1

        # Check if we need to move to the next ring.clear
        if self.current_ring == -1:
            self.current_ring += 1
        elif self.current_tract > sum(self.ring_nums[: self.current_ring + 1]):
            self.current_ring += 1

        return tract_and_ring


def get_poly_from_tract_id(skymap, tract_id, inner=False):
    """Get the ConvexPolygon for a tract by its ID.

    Parameters
    ----------
    skymap : lsst.skymap.SkyMap
        The LSST SkyMap object.
    tract_id : int
        The ID of the tract to retrieve.
    inner : bool, optional
        If True, return the inner polygon. If False, return the outer polygon.
        Default is False (outer polygon).

    Returns
    -------
    ConvexPolygon
        The polygon representing the tract's sky region.
    """
    from lsst.sphgeom import Box

    tract = skymap.generateTract(tract_id)
    res = tract.inner_sky_region if inner else tract.outer_sky_polygon
    res = box_to_convex_polygon(res) if isinstance(res, Box) else res
    return res


def polys_are_equiv(poly_a, poly_b, rtol=1e-12, atol=1e-14):
    """Check if two ConvexPolygons are equivalent within floating point tolerance.

    Parameters
    ----------
    poly_a, poly_b : sphgeom.ConvexPolygon
        The polygons to compare.
    rtol : float
        Relative tolerance for np.allclose.
    atol : float
        Absolute tolerance for np.allclose.

    Returns
    -------
    bool
        True if all vertices match within tolerance.
    """
    verts_a = poly_a.getVertices()
    verts_b = poly_b.getVertices()

    if len(verts_a) != len(verts_b):
        return False

    return np.allclose(verts_a, verts_b, rtol=rtol, atol=atol)


def quads_are_equiv(
    quad1: list[list[float]],
    quad2: list[list[float]],
    tol_arcsec: float = 1.0,
) -> bool:
    """Check if two quads are approximately equivalent in RA/Dec space.

    Parameters
    ----------
    quad1 : list of [RA, Dec]
        First set of four polygon vertices in degrees.
    quad2 : list of [RA, Dec]
        Second set of four polygon vertices in degrees.
    tol_arcsec : float, optional
        Allowed tolerance in arcseconds for all matching points.

    Returns
    -------
    bool
        True if all corresponding vertices are within the specified tolerance.
    """
    if len(quad1) != 4 or len(quad2) != 4:
        raise ValueError("Each quad must contain exactly 4 vertices.")

    # Convert to SkyCoord objects
    c1 = SkyCoord(ra=[ra for ra, dec in quad1] * u.deg, dec=[dec for ra, dec in quad1] * u.deg)
    c2 = SkyCoord(ra=[ra for ra, dec in quad2] * u.deg, dec=[dec for ra, dec in quad2] * u.deg)

    # Compute angular separation for each corresponding vertex
    separations = c1.separation(c2).arcsecond

    return np.all(separations < tol_arcsec)


def get_patch_poly_from_ids(skymap, tract_id, patch_id):
    """Get ConvexPolygon for a given patch."""
    patch_info = skymap.generateTract(tract_id).getPatchInfo(patch_id)
    return patch_info.inner_sky_polygon


def get_quad_from_tract_id(skymap, tract_id, inner=True) -> list[list[float]]:
    """Get the RA/Dec quad vertices for a tract by its ID.

    Parameters
    ----------
    skymap : lsst.skymap.SkyMap
        The LSST SkyMap object.
    tract_id : int
        The ID of the tract to retrieve.
    inner : bool, optional
        If True, return the inner quad. If False, return the outer quad.
        Default is True (inner).

    Returns
    -------
    list of [RA, Dec]
        The list of four polygon vertices.
    """
    from lsst.sphgeom import Box

    tract = skymap.generateTract(tract_id)
    res = tract.inner_sky_region if inner else tract.outer_sky_polygon
    res = box_to_convex_polygon(res) if isinstance(res, Box) else res
    return [unit_vector3d_to_radec(v) for v in res.getVertices()]


def get_quad_from_patch_id(skymap, tract_id: int, patch_id: int) -> list[list[float]]:
    """Get the RA/Dec quad vertices for a specific patch in a tract.

    Parameters
    ----------
    skymap : lsst.skymap.SkyMap
        The LSST SkyMap object.
    tract_id : int
        The ID of the tract containing the patch.
    patch_id : int
        The sequential ID of the patch (typically 0–99).

    Returns
    -------
    list of [RA, Dec]
        The list of four polygon vertices for the patch.
    """
    tract = skymap.generateTract(tract_id)
    patch_info = tract.getPatchInfo(patch_id)
    return [unit_vector3d_to_radec(v) for v in patch_info.inner_sky_polygon.getVertices()]
