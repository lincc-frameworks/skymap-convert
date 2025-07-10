import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from lsst.sphgeom import Box, ConvexPolygon

from .utils import box_to_convex_polygon, unit_vector3d_to_radec


def get_poly_from_tract_id(skymap, tract_id, inner=False) -> ConvexPolygon:
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
