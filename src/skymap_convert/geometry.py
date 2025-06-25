"""Geometry utilities for converting between different representations of sky regions."""

from lsst.sphgeom import ConvexPolygon, LonLat, UnitVector3d


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
    lonlat = LonLat(vec)
    ra = lonlat.getLon().asDegrees() % 360.0
    dec = lonlat.getLat().asDegrees()
    return ra, dec
