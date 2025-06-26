import numpy as np
from lsst.sphgeom import Box, ConvexPolygon

from .utils import box_to_convex_polygon


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
