import numpy as np
from lsst.sphgeom import Box, ConvexPolygon

from .geometry import box_to_convex_polygon


def get_poly_from_tract_id(skymap, tract_id, inner=False) -> ConvexPolygon:
    tract = skymap.generateTract(tract_id)
    if inner:
        res = tract.inner_sky_region
    else:
        res = tract.outer_sky_polygon
    if isinstance(res, Box):
        res = box_to_convex_polygon(res)
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
