import gzip
import math
import shutil
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path

import numpy as np
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
            The output file path or directory where the output file will be written.
        """
        if isinstance(output_path, str):
            output_path = Path(output_path)

        # If output_path is a file, get its parent directory
        if output_path.is_file():
            output_path = output_path.parent

        # Create the directory if it doesn't exist
        output_path.mkdir(parents=True, exist_ok=True)

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

        print(f"Wrote {len(out['tracts'])} tract polygons to {output_path}")


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
            print(f"Ring-optimized skymap written to {output_path}")


class ConvertedSkymapWriter(SkymapWriter):
    """Writer for Converted Skymaps using .npy + metadata.yaml."""

    def write(self, skymap, output_path: str | Path, skymap_name: str = "converted_skymap"):
        """Write tract and patch vertices to .npy files with metadata.

        Parameters
        ----------
        skymap : lsst.skymap.SkyMap
            The LSST SkyMap object to write
        output_path : str or Path
            Directory path to write output files into
        skymap_name : str, optional
            Name of the skymap (for metadata purposes)
        """
        output_path = Path(output_path)
        self._ensure_output_directory(output_path)

        n_tracts = len(skymap)  # TODO - is this accounting for poles?
        n_patches = 100  # fixed per tract # TODO - would be nice to make this dynamic, though we expect 100
        tract_array = np.zeros((n_tracts, 4, 2), dtype=np.float64)
        patch_array = np.zeros((n_tracts, n_patches, 4, 2), dtype=np.float64)

        for tract_data in skymap:
            tract_id = tract_data.getId()

            # Get tract vertices (inner region only)
            inner_sky_region = tract_data.inner_sky_region
            if isinstance(inner_sky_region, Box):
                inner_sky_region = box_to_convex_polygon(inner_sky_region)
            tract_verts = [unit_vector3d_to_radec(vert) for vert in inner_sky_region.getVertices()]
            tract_array[tract_id] = tract_verts

            # Get patch vertices
            for patch_id in range(n_patches):
                patch_info = tract_data.getPatchInfo(patch_id)
                patch_verts = [
                    unit_vector3d_to_radec(vert) for vert in patch_info.inner_sky_polygon.getVertices()
                ]
                patch_array[tract_id, patch_id] = patch_verts

        # Save tracts (uncompressed)
        np.save(output_path / "tracts.npy", tract_array)

        # Save patches (temporary .npy before compression)
        patches_path = output_path / "patches.npy"
        patches_gz_path = output_path / "patches.npy.gz"

        np.save(patches_path, patch_array)

        # Compress to patches.npy.gz
        with open(patches_path, "rb") as f_in, gzip.open(patches_gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Remove uncompressed version
        patches_path.unlink()

        # Save metadata
        metadata = {
            "name": skymap_name,
            "generated": datetime.utcnow().isoformat() + "Z",
            "n_tracts": n_tracts,
            "n_patches_per_tract": n_patches,
            "format_version": 1,
        }
        with open(output_path / "metadata.yaml", "w") as f:
            yaml.dump(metadata, f)

        print(f"Wrote Converted Skymap to: {output_path}")
