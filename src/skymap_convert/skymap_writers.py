import math
from abc import ABC, abstractmethod
from pathlib import Path

import yaml
from lsst.sphgeom import Box

from .utils import box_to_convex_polygon, radians_to_degrees, unit_vector3d_to_radec


class SkymapWriter(ABC):
    """Abstract base class for writing skymaps to files."""

    def _ensure_output_directory(self, output_path: Path):
        """Ensure the output directory exists.

        Parameters
        ----------
        output_path : Path
            The output file path
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def write(self, skymap, output_path: str | Path, **kwargs):
        """Write the skymap to file.

        Parameters
        ----------
        skymap : lsst.skymap.SkyMap
            The LSST SkyMap object to write
        output_path : str or Path
            Destination path for the output file
        **kwargs
            Additional keyword arguments specific to each writer
        """
        pass


class FullVertexWriter(SkymapWriter):
    """Writer for full vertex format skymaps."""

    def write(self, skymap, output_path: str | Path, inner=True, write_patches=False):
        """Write tract (and optionally patch) polygons to YAML using RA/Dec coordinates.

        Parameters
        ----------
        skymap : lsst.skymap.SkyMap
            The LSST SkyMap object.
        output_path : str or Path
            Destination path for the output file
        inner : bool, optional
            If True, write inner polygons. If False, write outer polygons. Default is True.
        write_patches : bool, optional
            If True, include patch polygons for each tract. Default is False.
        """
        output_path = Path(output_path)
        out = {"tracts": {}}

        for tract in skymap:
            tract_id = tract.getId()
            if inner:
                poly = tract.inner_sky_region
                if isinstance(poly, Box):
                    poly = box_to_convex_polygon(poly)
            else:
                poly = tract.outer_sky_polygon

            ra_dec_vertices = [
                [radec[0], radec[1]] for radec in map(unit_vector3d_to_radec, poly.getVertices())
            ]
            out["tracts"][tract_id] = {"polygon": ra_dec_vertices}

            # Save patch vectors, too.
            # TODO, but don't; we'll use a different implementation.

        # Ensure output directory exists
        self._ensure_output_directory(output_path)

        # Write to YAML file.
        with open(output_path, "w") as f:
            yaml.dump(out, f, sort_keys=False)

        print(f"✅ Wrote {len(out['tracts'])} tract polygons to {output_path}")


class RingOptimizedWriter(SkymapWriter):
    """Writer for ring-optimized format skymaps."""

    def write(
        self,
        skymap,
        output_path: str | Path,
        inner=True,
        patches=False,
        skymap_name=None,
    ):
        """Write a ring-optimized skymap to YAML format.

        Parameters
        ----------
        skymap : lsst.skymap.SkyMap
            The LSST SkyMap object.
        output_path : str or Path
            Destination path for the output file
        inner : bool, optional
            If True, include inner polygons in the output. Default is True.
        patches : bool, optional
            If True, include patch polygons in the output. Default is False.
        skymap_name : str, optional
            Name of the skymap, used in the metadata. If None, defaults to "ring_optimized_skymap".
        """
        output_path = Path(output_path)
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
        # todo : patch metadata, like how many per tract

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

        else:
            raise NotImplementedError(
                "Outer polygons are not yet implemented for ring-optimized skymaps. "
                "Please use Full Vertex skymap for outer polygons instead."
            )

        # Ensure output directory exists
        self._ensure_output_directory(output_path)

        # Record the output.
        with open(output_path, "w") as f:
            yaml.dump(out, f, sort_keys=False)
            print(f"✅ Ring-optimized skymap written to {output_path}")
