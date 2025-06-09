"""Utility functions for reading and writing skymap polygons in RA/Dec coordinates."""

import yaml
import pickle
from pathlib import Path
from lsst.sphgeom import UnitVector3d, ConvexPolygon
from .geometry import box_to_convex_polygon, unit_vector3d_to_radec


def load_pickle_skymap(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def write_polygons_ra_dec(skymap, output_path, inner=True, write_patches=False):
    """Write tract (and optionally patch) polygons to YAML using RA/Dec coordinates.

    Parameters
    ----------
    skymap : `lsst.skymap.SkyMap`
        The LSST SkyMap object.
    output_path : str or Path
        Destination path for the output YAML file.
    inner : bool, optional
        If True, write inner polygons. If False, write outer polygons. Default is True.
    write_patches : bool, optional
        If True, include patch polygons for each tract. Default is False.
    """
    import lsst.geom as geom
    from lsst.sphgeom import LonLat, Box

    out = {"tracts": {}}

    for tract in skymap:
        tract_id = tract.getId()
        if inner:
            poly = tract.inner_sky_region
            if isinstance(poly, Box):
                poly = box_to_convex_polygon(poly)
        else:
            poly = tract.outer_sky_polygon

        ra_dec_vertices = [[radec[0], radec[1]] for radec in map(unit_vector3d_to_radec, poly.getVertices())]
        out["tracts"][tract_id] = {"polygon": ra_dec_vertices}

        # Save patch vectors, too.
        # NB: The theory is here, and this technically works, but will be much better to do this
        # after implementing the ring-based optimization.
        # if write_patches:
        #     patch_polys = []
        #     for patch in tract:
        #         corners = patch.getOuterBBox().getCorners()
        #         wcs = tract.getWcs()
        #         patch_ra_dec = []
        #         for corner in corners:
        #             sky_pt = wcs.pixelToSky(geom.Point2D(corner))
        #             ra = sky_pt.getLongitude().asDegrees() % 360.0
        #             dec = sky_pt.getLatitude().asDegrees()
        #             patch_ra_dec.append([ra, dec])
        #         patch_polys.append(patch_ra_dec)
        #     out["tracts"][tract_id]["patches"] = patch_polys

    with open(output_path, "w") as f:
        yaml.dump(out, f, sort_keys=False)


# def load_polygons_ra_dec(yaml_path):
#     """Load tract polygons from a YAML file written with RA/Dec coordinates.

#     Parameters
#     ----------
#     yaml_path : str or Path
#         Path to the YAML file created by `write_polygons_ra_dec`.

#     Returns
#     -------
#     dict
#         Dictionary mapping tract ID (int) to `sphgeom.ConvexPolygon`.
#     """
#     from lsst.sphgeom import LonLat, UnitVector3d, ConvexPolygon

#     with open(yaml_path, "r") as f:
#         data = yaml.safe_load(f)

#     poly_dict = {}
#     for tract_id_str, record in data["tracts"].items():
#         tract_id = int(tract_id_str)
#         vertices = record["polygon"]

#         unit_vecs = [UnitVector3d.fromLonLat(*pair) for pair in vertices]
#         unique_vecs = {tuple(round(x, 12) for x in v) for v in unit_vecs}
#         if len(unique_vecs) < 3:
#             print(f"⚠️ Storing `None` for degenerate tract {tract_id}")
#             poly_dict[tract_id] = None
#             continue
#         poly_dict[tract_id] = ConvexPolygon(unit_vecs)

#     return poly_dict

from lsst.sphgeom import LonLat, UnitVector3d, ConvexPolygon


def load_polygons_ra_dec(yaml_path):
    """Load tract polygons from a YAML file written with RA/Dec coordinates.

    Parameters
    ----------
    yaml_path : str or Path
        Path to the YAML file created by `write_polygons_ra_dec`.

    Returns
    -------
    dict
        Dictionary mapping tract ID (int) to `sphgeom.ConvexPolygon` or None if degenerate.
    """
    import yaml

    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)

    poly_dict = {}

    for tract_id_str, content in data["tracts"].items():
        tract_id = int(tract_id_str)
        ra_dec_vertices = content["polygon"]

        unit_vecs = [UnitVector3d(LonLat.fromDegrees(ra % 360.0, dec)) for ra, dec in ra_dec_vertices]

        # Round for precision-safe uniqueness check
        unique_vecs = {tuple(round(coord, 12) for coord in vec) for vec in unit_vecs}

        if len(unique_vecs) < 3:
            print(f"⚠️ Storing `None` for degenerate tract {tract_id}")
            poly_dict[tract_id] = None
            continue

        poly_dict[tract_id] = ConvexPolygon(unit_vecs)

    print(f"✅ Loaded {len(poly_dict)} tract polygons from {yaml_path}")
    return poly_dict
