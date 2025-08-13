import gzip
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

from .utils import box_to_convex_polygon, unit_vector3d_to_radec


class ConvertedSkymapWriter:
    """Writer for Converted Skymaps using .npy + metadata.yaml format.

    This class converts LSST SkyMap objects into a format optimized for fast loading
    and analysis. It extracts tract and patch polygon vertices and saves them as
    NumPy arrays with accompanying YAML metadata.

    Notes
    -----
    This class requires the `lsst.sphgeom` package.

    Examples
    --------
    >>> writer = ConvertedSkymapWriter()
    >>> writer.write(skymap, "output_dir", "my_skymap")
    """

    def write(self, skymap, output_path: str | Path, skymap_name: str = "converted_skymap"):
        """Write tract and patch vertices to .npy files with metadata.

        Extracts polygon vertices from an LSST SkyMap and saves them in an optimized
        format consisting of:
        - tracts.npy: Uncompressed array of tract boundary vertices
        - patches.npy.gz: Compressed array of patch boundary vertices
        - metadata.yaml: Skymap information and array dimensions

        Parameters
        ----------
        skymap : lsst.skymap.SkyMap
            The LSST SkyMap object to convert and write.
        output_path : str or Path
            Directory path where output files will be written.
        skymap_name : str, optional
            Descriptive name for the skymap (stored in metadata), by default "converted_skymap".

        Notes
        -----
        The method assumes 100 patches per tract, which is standard for LSST skymaps.
        All polygon vertices are stored as [RA, Dec] coordinates in degrees.

        Examples
        --------
        >>> writer = ConvertedSkymapWriter()
        >>> writer.write(my_skymap, "/path/to/output", "DC2_skymap")
        """
        from lsst.sphgeom import Box

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

    def _ensure_output_directory(self, output_path: Path):
        """Ensure the output directory exists, creating it if necessary.

        Parameters
        ----------
        output_path : Path
            The directory path where output files will be written.

        Notes
        -----
        Creates parent directories as needed. If output_path points to an existing
        file, uses the file's parent directory instead.
        """
        if isinstance(output_path, str):
            output_path = Path(output_path)

        # If output_path is a file, get its parent directory
        if output_path.is_file():
            output_path = output_path.parent

        # Create the directory if it doesn't exist
        output_path.mkdir(parents=True, exist_ok=True)
