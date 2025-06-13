"""Utility functions for reading and writing skymap polygons in RA/Dec coordinates."""

import math
import pickle
from pathlib import Path

import yaml
from lsst.sphgeom import ConvexPolygon, LonLat, UnitVector3d

from .geometry import box_to_convex_polygon, unit_vector3d_to_radec
from .utils import radians_to_degrees


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
    from lsst.sphgeom import Box, LonLat

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

    # If output directory does not exist, create it
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write to YAML file.
    with open(output_path, "w") as f:
        yaml.dump(out, f, sort_keys=False)

    print(f"✅ Wrote {len(out['tracts'])} tract polygons to {output_path}")


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


def write_ring_optimized_skymap(skymap, output_path, inner=True, patches=False, skymap_name=None):
    """Write a ring-optimized skymap to YAML format.

    Parameters
    ----------
    skymap : `lsst.skymap.SkyMap`
        The LSST SkyMap object.
    output_path : str or Path
        Destination path for the output YAML file.
    inner : bool, optional
        If True, include inner polygons in the output. Default is True.
    patches : bool, optional
        If True, include patch polygons in the output. Default is False.
    skymap_name : str, optional
        Name of the skymap, used in the metadata. If None, defaults to "ring_optimized_skymap".
    """
    ring_nums = skymap._ringNums
    ring_size = math.pi / (len(ring_nums) + 1)
    total_tracts = sum(ring_nums) + 2  # +2 for the poles

    # Initialize the output structure.
    out = {"metadata": {}, "poles": [], "rings": []}

    # Add the metadata.
    if skymap_name is None:
        skymap_name = "ring_optimized_skymap"
    out["metadata"]["skymap_name"] = skymap_name
    out["metadata"]["inner_polygons"] = inner
    out["metadata"]["ra_start"] = skymap.config.raStart

    # todo : read through how many patches in each tract,
    # or is this into at a higher level? then include

    # Handle the poles first.
    first_ring_middle_dec = ring_size * (1) - 0.5 * math.pi
    first_ring_lower_dec = radians_to_degrees(first_ring_middle_dec - 0.5 * ring_size)
    south_pole = {
        "tract_id": 0,
        "ring": -1,
        "dec_bounds": [-90.0, first_ring_lower_dec],
        "ra_bounds": [0.0, 360.0],
        "ra_interval": 360.0,
    }
    out["poles"].append(south_pole)
    last_ring_middle_dec = ring_size * len(ring_nums) - 0.5 * math.pi
    last_ring_upper_dec = radians_to_degrees(last_ring_middle_dec + 0.5 * ring_size)
    north_pole = {
        "tract_id": total_tracts - 1,
        "ring": len(ring_nums),
        "dec_bounds": [last_ring_upper_dec, 90.0],
        "ra_bounds": [0.0, 360.0],
        "ra_interval": 360.0,
    }
    out["poles"].append(north_pole)

    # Now the rings.
    if inner:
        # Add the rings.
        tract_counter = 1  # tract 0 is pole. this var is temp, for early breaking.
        for ring, num_tracts in enumerate(ring_nums):
            # Get the declination bounds for the ring.
            dec = ring_size * (ring + 1) - 0.5 * math.pi
            start_dec = radians_to_degrees(dec - 0.5 * ring_size)
            stop_dec = radians_to_degrees(dec + 0.5 * ring_size)

            # Get the RA interval for the ring.
            ra_interval = 360.0 / num_tracts

            # Write and add the ring entry.
            ring_entry = {
                "ring": ring,
                "num_tracts": num_tracts,
                "dec_bounds": [round(start_dec, 10), round(stop_dec, 10)],
                "ra_interval": round(ra_interval, 10),
            }
            out["rings"].append(ring_entry)

            # Iterate the tract counter.
            tract_counter += num_tracts

            # Temporary: just break after 10 rings for testing.
            # if ring >= 10:
            #     break

    else:
        print("oh no!")

    # If directory does not exist, create it.
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # Record the output.
    with open(output_path, "w") as f:
        yaml.dump(out, f, sort_keys=False)
        print(f"✅ Ring-optimized skymap written to {output_path}")


class RingOptimizedSkymapReader:
    """A reader for ring-optimized skymaps written in YAML format.

    This class reads the YAML file and provides access to the metadata, rings, tracts, and poles.
    """

    def __init__(self, yaml_path):
        with open(yaml_path, "r") as f:
            self.data = yaml.safe_load(f)

        self.metadata = self.data["metadata"]
        self.skymap_name = self.metadata.get("skymap_name", "ring_optimized_skymap")
        self.inner_polygons = self.metadata.get("inner_polygons", True)
        self.rings = self.data["rings"]
        self.poles = self.data["poles"]
        self.ra_start = self.metadata.get("ra_start", 0.0)
        self.ring_nums = [ring["num_tracts"] for ring in self.rings]
        self.total_tracts = sum(self.ring_nums) + len(self.poles)

    def __str__(self):
        return (
            f"RingOptimizedSkymapReader("
            f"skymap_name={self.skymap_name}, "
            f"rings={len(self.rings)}, "
            f"poles={len(self.poles)}, "
            f"inner_polygons={self.inner_polygons}"
            f")"
        )

    def get(self, tract_id):
        """A wrapper for `get_tract` to allow dictionary-like access.

        Parameters
        ----------
        tract_id : int
            The index of the tract to retrieve.

        Returns
        -------
        dict or None
            The tract data, or None if the index is out of bounds.
        """
        return self.get_tract(tract_id)

    def get_ring(self, ring_index):
        """Get the ring data for a given ring index.

        Parameters
        ----------
        ring_index : int
            The index of the ring to retrieve.

        Returns
        -------
        dict
            The ring data, or None if the index is out of bounds.
        """
        if 0 <= ring_index < len(self.rings):
            return self.rings[ring_index]
        else:
            return None

    def get_ring_from_tract_id(self, tract_index):
        """Get the ring data for a given tract index.

        Parameters
        ----------
        tract_index : int
            The index of the tract to retrieve the ring for.

        Returns
        -------
        dict or None
            The ring data, or None if the tract index is out of bounds.
        """
        if 0 <= tract_index < self.total_tracts:
            # Find the ring for the given tract index.
            for ring in self.rings:
                if tract_index < sum(self.ring_nums[: ring["ring"] + 1]):
                    return ring
            return None
        else:
            return None

    def get_pole(self, which_pole):
        """Get the pole data for the specified pole.

        Parameters
        ----------
        which_pole : str
            Either "north" or "south".

        Returns
        -------
        dict or None
            The pole data, or None if the pole is not defined.
        """
        if which_pole.lower() == "north":
            return self.poles[1] if len(self.poles) > 1 else None
        elif which_pole.lower() == "south":
            return self.poles[0] if self.poles else None
        else:
            raise ValueError("which_pole must be 'north' or 'south'.")

    def _construct_quad(self, dec_min, dec_max, ra_start, ra_end, pole=False):
        """Construct a quadrilateral polygon from the given bounds.

        Quadrilateral is defined by four vertices in RA/Dec coordinates (degrees). Vertices begin at
        the lower left corner and go counter-clockwise.

        Note that, with RA wrapping, the quadrilateral may cross the RA=0 line, so the RA values may
        not be in increasing order.

        Parameters
        ----------
        dec_min : float
            Minimum declination in degrees.
        dec_max : float
            Maximum declination in degrees.
        ra_start : float
            Starting right ascension in degrees.
        ra_end : float
            Ending right ascension in degrees.
        pole : bool, optional
            If True, the quadrilateral is constructed for a pole (south or north). Default is False.
            If True, the RA values are not wrapped to [0, 360) range (to allow for an RA = 360).
            todo this may change

        Returns
        -------
        list of list of float
            A list of [RA, Dec] representing the vertices of the quadrilateral.
        """
        if pole:
            return [
                [ra_start, dec_min],  # Lower left
                [ra_end, dec_min],  # Lower right
                [ra_end, dec_max],  # Upper right
                [ra_start, dec_max],  # Upper left
            ]
        return [
            # (ra_start % 360.0, dec_min),  # Lower left
            # (ra_end % 360.0, dec_min),  # Lower right
            # (ra_end % 360.0, dec_max),  # Upper right
            # (ra_start % 360.0, dec_max),  # Upper left
            [ra_start, dec_min],  # Lower left
            [ra_end, dec_min],  # Lower right
            [ra_end, dec_max],  # Upper right
            [ra_start, dec_max],  # Upper left
        ]

    def _construct_tract_data(self, tract_id, ring=None):
        """Construct the tract data dictionary.

        Parameters
        ----------
        tract_id : int
            The ID of the tract.
        ring :  dict
            The ring data for the tract, if tract is not a pole.

        Returns
        -------
        dict
            A dictionary containing the tract data.
        """
        # Out of bounds check.
        if tract_id < 0 or tract_id >= self.total_tracts:
            raise ValueError(f"Tract ID {tract_id} is out of bounds for this skymap.")

        # Handle the poles. (todo) would be nice to have this as a separate method.
        elif tract_id == 0:
            if self.poles:
                dec_min, dec_max = self.poles[0]["dec_bounds"]
                ra_start, ra_end = self.poles[0]["ra_bounds"][0], self.poles[0]["ra_bounds"][1]
                ra_start = ra_start % 360.0
                ra_end = ra_end % 360.0
                quad = self._construct_quad(dec_min, dec_max, ra_start, ra_end, pole=False)
                return {
                    "tract_id": tract_id,
                    "ring": -1,
                    "dec_bounds": [dec_min, dec_max],
                    "ra_bounds": [ra_start, ra_end],
                    "quad": quad,
                }
            else:
                raise ValueError("No south pole defined in the skymap.")
        elif tract_id == self.total_tracts - 1:
            if self.poles and len(self.poles) > 1:
                dec_min, dec_max = self.poles[1]["dec_bounds"]
                ra_start, ra_end = self.poles[1]["ra_bounds"][0], self.poles[1]["ra_bounds"][1]
                ra_start = ra_start % 360.0
                ra_end = ra_end % 360.0
                quad = self._construct_quad(dec_min, dec_max, ra_start, ra_end, pole=True)
                return {
                    "tract_id": tract_id,
                    "ring": len(self.rings),
                    "dec_bounds": [dec_min, dec_max],
                    "ra_bounds": [ra_start, ra_end],
                    "quad": quad,
                }
            else:
                raise ValueError("No north pole defined in the skymap.")

        # Normal tracts in rings.
        dec_min, dec_max = ring["dec_bounds"]
        ra_interval = ring["ra_interval"]
        ra_start = self.ra_start + (tract_id - 1) * ra_interval - ra_interval * 0.5
        ra_end = ra_start + ra_interval
        ra_start = ra_start % 360.0
        ra_end = ra_end % 360.0
        quad = self._construct_quad(dec_min, dec_max, ra_start, ra_end)
        return {
            "tract_id": tract_id,
            "ring": ring,
            "dec_bounds": [dec_min, dec_max],
            "ra_bounds": [ra_start, ra_end],
            "quad": quad,
        }

    def get_tract(self, tract_id):
        """Get the tract data for the skymap.

        Parameters
        ----------
        tract_id : int
            The ID of the tract to retrieve.

        Returns
        -------
        dict
            The tract data.

        Raises
        ------
        ValueError
            If the tract ID is out of bounds or not found in the skymap.
        """
        # Handle the poles.
        if tract_id == 0 or tract_id == self.total_tracts - 1:
            return self._construct_tract_data(tract_id)

        # Handle normal tracts in rings.
        # tract_id -= 1  # Adjust for the south pole being tract 0 and being unrepresented in self.rings.
        # print(f"Getting tract {tract_id}...")
        for ring in self.rings:
            tracts_in_ring = ring["num_tracts"]
            if tract_id <= tracts_in_ring:  # tract is in this ring
                # print("Found tract in ring:", ring["ring"])
                return self._construct_tract_data(tract_id, ring)
            else:
                tract_id -= tracts_in_ring

        return ValueError(f"Tract ID {tract_id} not found in the skymap.")
