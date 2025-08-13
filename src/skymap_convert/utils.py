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
    lsst.sphgeom.ConvexPolygon
        The equivalent polygon with four vertices.

    Raises
    ------
    ValueError
        If the box is empty.

    Notes
    -----
    The conversion assumes the box does not wrap around the RA=0/360 boundary.
    This may need adjustment for boxes that span this boundary.

    TODO : Handle RA wrap-around cases more robustly.
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
    """Convert a UnitVector3d to RA/Dec coordinates in degrees.

    Parameters
    ----------
    vec : lsst.sphgeom.UnitVector3d
        A 3D unit vector representing a point on the sphere.

    Returns
    -------
    tuple of float
        (RA, Dec) coordinates where RA is in range [0, 360) degrees
        and Dec is in range [-90, 90] degrees.
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

    Raises
    ------
    FileNotFoundError
        If the pickle file does not exist.
    pickle.UnpicklingError
        If the file cannot be unpickled or does not contain a valid SkyMap.
    """
    with open(path, "rb") as f:
        return pickle.load(f)


def radians_to_degrees(radians):
    """Convert radians to degrees.

    Parameters
    ----------
    radians : float
        Angle in radians.

    Returns
    -------
    float
        Angle in degrees.
    """
    return radians * (180.0 / math.pi)


class IterateTractAndRing:
    """An iterator to traverse through tract IDs and ring IDs in a skymap.

    This iterator yields tuples of (tract_index, ring_index) for each tract in the skymap.
    The first tract is the south pole (ring -1), and the last tract is the north pole.

    Parameters
    ----------
    ring_nums : list of int
        A list where each element represents the number of tracts in each ring.
        For example, [5, 10, 15] represents 5 tracts in ring 0, 10 in ring 1,
        and 15 in ring 2.
    add_poles : bool, optional
        If True, include the south pole (tract 0) and north pole (last tract)
        in the iteration. If False, only iterate through the rings (default: True).

    Examples
    --------
    >>> iterator = IterateTractAndRing([5, 10, 15], add_poles=True)
    >>> for tract_id, ring_id in iterator:
    ...     print(f"Tract {tract_id} is in ring {ring_id}")
    Tract 0 is in ring -1
    Tract 1 is in ring 0
    ...

    TODO : this is more or less deprecated, I belive. Consider removing.
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
        If True, return the inner polygon. If False, return the outer polygon
        (default: False).

    Returns
    -------
    lsst.sphgeom.ConvexPolygon
        The polygon representing the tract's sky region.

    Notes
    -----
    Automatically handles conversion from Box to ConvexPolygon if needed.
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
    poly_a, poly_b : lsst.sphgeom.ConvexPolygon
        The polygons to compare.
    rtol : float, optional
        Relative tolerance for numpy.allclose (default: 1e-12).
    atol : float, optional
        Absolute tolerance for numpy.allclose (default: 1e-14).

    Returns
    -------
    bool
        True if all vertices match within tolerance, False otherwise.
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
    quad1 : list of list of float
        First set of four polygon vertices as [[RA, Dec], ...] in degrees.
    quad2 : list of list of float
        Second set of four polygon vertices as [[RA, Dec], ...] in degrees.
    tol_arcsec : float, optional
        Allowed tolerance in arcseconds for all matching points (default: 1.0).

    Returns
    -------
    bool
        True if all corresponding vertices are within the specified tolerance.

    Raises
    ------
    ValueError
        If either quad does not contain exactly 4 vertices.
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
    """Get ConvexPolygon for a given patch.

    Parameters
    ----------
    skymap : lsst.skymap.SkyMap
        The LSST SkyMap object.
    tract_id : int
        The ID of the tract containing the patch.
    patch_id : int
        The ID of the patch within the tract.

    Returns
    -------
    lsst.sphgeom.ConvexPolygon
        The polygon representing the patch's inner sky region.
    """
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
        If True, return the inner quad. If False, return the outer quad
        (default: True).

    Returns
    -------
    list of list of float
        List of four polygon vertices as [[RA, Dec], ...] in degrees.

    Notes
    -----
    Automatically handles conversion from Box to ConvexPolygon if needed.
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
    list of list of float
        List of four polygon vertices as [[RA, Dec], ...] in degrees.
    """
    tract = skymap.generateTract(tract_id)
    patch_info = tract.getPatchInfo(patch_id)
    return [unit_vector3d_to_radec(v) for v in patch_info.inner_sky_polygon.getVertices()]
