import numpy as np
from lsst.sphgeom import Box, ConvexPolygon, LonLat, UnitVector3d

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


def box_to_convex_polygon(box: Box) -> ConvexPolygon:
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

    bottom_left = LonLat.fromRadians(lon_min, lat_min)
    bottom_right = LonLat.fromRadians(lon_max, lat_min)
    top_right = LonLat.fromRadians(lon_max, lat_max)
    top_left = LonLat.fromRadians(lon_min, lat_max)

    # Convert corners to UnitVector3d
    vertices = [
        UnitVector3d(bottom_left),
        UnitVector3d(bottom_right),
        UnitVector3d(top_right),
        UnitVector3d(top_left),
    ]

    # Create and return the ConvexPolygon
    return ConvexPolygon(vertices)
