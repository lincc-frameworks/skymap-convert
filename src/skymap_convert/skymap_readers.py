import gzip
import shutil
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import yaml
from lsst.sphgeom import ConvexPolygon, LonLat, UnitVector3d

from .plotting import plot_patches
from .tract_data import TractData


class SkymapReader(ABC):
    """Abstract base class for reading skymaps from files."""

    def __init__(self, file_path: str | Path):
        """Initialize the reader with file path.

        Parameters
        ----------
        file_path : str or Path
            Path to the skymap file to read
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"Skymap file not found: {file_path}")

    def help(self):
        """Display available public methods and their descriptions.

        This method introspects the current reader instance to show all
        available public methods along with brief descriptions extracted
        from their docstrings.
        """
        import inspect

        class_name = self.__class__.__name__
        print(f"{class_name} - Available Methods")
        print("=" * (len(class_name) + 21))
        print()

        # Get all methods that don't start with underscore (public methods)
        methods = []
        for name, method in inspect.getmembers(self, predicate=inspect.ismethod):
            if not name.startswith("_"):
                methods.append((name, method))

        # Sort methods alphabetically
        methods.sort(key=lambda x: x[0])

        for method_name, method in methods:
            # Get the first line of the docstring as brief description
            doc = inspect.getdoc(method)
            if doc:
                # Take first sentence or first line as brief description
                brief = doc.split(".")[0] + "." if "." in doc else doc.split("\n")[0]
                # Limit length to keep it concise
                if len(brief) > 80:
                    brief = brief[:77] + "..."
            else:
                brief = "No description available"

            print(f"  {method_name:20} - {brief}")

        print()
        print("For detailed information about any method, use: help(reader.method_name)")

    @abstractmethod
    def get_tract(self, tract_id: int) -> TractData:
        """Get tract data for a given tract ID (legacy; may be deprecated).

        Parameters
        ----------
        tract_id : int
            The tract ID to retrieve

        Returns
        -------
        TractData
            The tract data
        """
        pass

    # @abstractmethod
    # def get_all_tracts(self) -> dict:
    #     """Get all tract data.

    #     Returns
    #     -------
    #     dict
    #         Dictionary mapping tract ID to TractData
    #     """
    #     pass


class ConvertedSkymapReader(SkymapReader):
    """Reader for Converted Skymaps written as .npy files with metadata.

    TODO : make a list of the attributes people might want
    - n_tracts : int
    - n_patches_per_tract : int
    - metadata : dict
    - tracts : np.ndarray
    - patches : np.ndarray
    """

    def __init__(self, file_path: str | Path = None, safe_loading: bool = False, preset: str = None):
        """Initialize the reader and load the .npy + metadata files.

        Parameters
        ----------
        file_path : str or Path, optional
            Path to the directory containing tracts.npy, patches.npy, and metadata.yaml.
            If None, must specify a preset.
        safe_loading : bool, optional
            If True, raise an exception on degenerate tract or patch polygons.
        preset : str, optional
            Name of a built-in skymap preset to load. If specified, file_path is ignored.
            Available presets can be listed with skymap_convert.presets.list_available_presets().
        """
        # Handle preset parameter
        if preset is not None:
            from .presets import get_preset_path

            file_path = get_preset_path(preset)

        elif file_path is None:
            raise ValueError("Either file_path or preset must be specified")

        super().__init__(file_path)
        self.safe_loading = safe_loading

        # Check if the path is a pickle file (eg, a pre-converted skymap)
        if self.file_path.is_file() and self.file_path.suffix in [".pickle", ".pkl"]:
            raise ValueError(
                f"ConvertedSkymapReader expects a directory containing converted skymap files, "
                f"but received a pickle file: {self.file_path}. "
                f"To read pickle files, first convert them using ConvertedSkymapWriter or "
                f"load them directly using load_pickle_skymap() from skymap_convert.utils."
            )

        self.metadata_path = self.file_path / "metadata.yaml"
        self.tracts_path = self.file_path / "tracts.npy"
        self.patches_path = self._decompress_patches_gz()
        self.patches = np.load(self.patches_path, mmap_mode="r")

        # Load metadata
        with open(self.metadata_path, "r") as f:
            self.metadata = yaml.safe_load(f)

        self.n_tracts = self.metadata["n_tracts"]
        self.n_patches_per_tract = self.metadata["n_patches_per_tract"]

        # Memory-map arrays
        self.tracts = np.load(self.tracts_path, mmap_mode="r")
        self.patches = np.load(self.patches_path, mmap_mode="r")

        # TODO could be nice to check if the metadata matches the arrays here.

    def _decompress_patches_gz(self) -> Path:
        """Decompress patches.npy.gz to a temp file if not already done."""
        patches_gz_path = self.file_path / "patches.npy.gz"

        # If already decompressed during this session, just return it.
        if hasattr(self, "_tmp_patches_path") and self._tmp_patches_path:
            return self._tmp_patches_path

        # Decompress to temp file
        with (
            gzip.open(patches_gz_path, "rb") as f_in,
            tempfile.NamedTemporaryFile(suffix=".npy", delete=False) as tmp_file,
        ):
            shutil.copyfileobj(f_in, tmp_file)
            self._tmp_patches_path = Path(tmp_file.name)

        return self._tmp_patches_path

    def cleanup(self):
        """Delete the temporary decompressed patches file, if it exists."""
        if hasattr(self, "_tmp_patches_path") and self._tmp_patches_path:
            try:
                self._tmp_patches_path.unlink()
            except Exception as e:
                print(f"Warning: Could not delete temp file {self._tmp_patches_path}: {e}")
            self._tmp_patches_path = None

    def _verify_nondegeneracy(self, vertices: list[list[float]]):
        """Verify that the polygon is non-degenerate by checking area.

        Parameters
        ----------
        vertices : list of [RA, Dec]
            The 4 polygon vertices (in degrees)

        Raises
        ------
        ValueError
            If polygon appears degenerate (zero-area)
        """
        if len(vertices) < 3:
            raise ValueError("Polygon has fewer than 3 vertices")

        # Use the first three vertices to compute area of the triangle they form
        a = np.array(vertices[0])
        b = np.array(vertices[1])
        c = np.array(vertices[2])

        # Vector AB and AC
        ab = b - a
        ac = c - a

        # Compute 2D cross product (scalar) to get area of parallelogram, then halve
        cross = ab[0] * ac[1] - ab[1] * ac[0]
        area = 0.5 * abs(cross)

        if area < 1e-6:  # Rough threshold for degeneracy at arcsecond scale
            raise ValueError("Degenerate polygon: near-zero area detected")

    def get_pole_tract_ids(self) -> list[int]:
        """Return the tract IDs that correspond to the poles.

        Returns
        -------
        list of int
            List of tract IDs for the poles
        """
        if not hasattr(self, "metadata") or "n_tracts" not in self.metadata:
            raise ValueError("Metadata is not loaded or does not contain 'n_tracts'")
        return [0, self.metadata["n_tracts"] - 1]

    def get_tract_vertices(self, tract_id: int) -> list[list[float]]:
        """Return the outer RA/Dec vertices of the specified tract.

        Parameters
        ----------
        tract_id : int
            ID of the tract to retrieve

        Returns
        -------
        list of [RA, Dec]
            List of four polygon vertices
        """
        if not (0 <= tract_id < self.n_tracts):
            raise IndexError(f"Tract ID {tract_id} is out of bounds")

        verts = self.tracts[tract_id].tolist()

        if self.safe_loading:
            self._verify_nondegeneracy(verts)

        return verts

    def get_patch_vertices(self, tract_id: int, patch_id: int) -> list[list[float]]:
        """Return the RA/Dec vertices of the specified patch.

        Parameters
        ----------
        tract_id : int
            ID of the tract containing the patch
        patch_id : int
            ID of the patch within the tract

        Returns
        -------
        list of [RA, Dec]
            List of four polygon vertices
        """
        if not (0 <= tract_id < self.n_tracts):
            raise IndexError(f"Tract ID {tract_id} is out of bounds")
        if not (0 <= patch_id < self.n_patches_per_tract):
            raise IndexError(f"Patch ID {patch_id} is out of bounds for tract {tract_id}")

        verts = self.patches[tract_id, patch_id].tolist()

        if self.safe_loading:
            self._verify_nondegeneracy(verts)

        return verts

    def get_tract(self, tract_id: int) -> TractData:
        """Return tract data for a given ID (slated for deprecation).

        Parameters
        ----------
        tract_id : int
            ID of the tract

        Returns
        -------
        TractData
            TractData object with quad and bounds

        Notes
        -----
        This method may be deprecated in future releases in favor of direct `get_tract_vertices`.
        """
        quad = self.get_tract_vertices(tract_id)

        ras = [ra for ra, _ in quad]
        decs = [dec for _, dec in quad]
        dec_bounds = (min(decs), max(decs))

        ra_min, ra_max = min(ras), max(ras)
        if ra_max - ra_min > 180:  # RA wraparound
            sorted_ras = sorted(ras)
            max_gap = 0
            gap_start = 0
            for i in range(len(sorted_ras)):
                gap = (sorted_ras[(i + 1) % 4] - sorted_ras[i]) % 360
                if gap > max_gap:
                    max_gap = gap
                    gap_start = sorted_ras[(i + 1) % 4]
            ra_bounds = (gap_start, (gap_start - max_gap) % 360)
        else:
            ra_bounds = (ra_min, ra_max)

        return TractData(
            tract_id=tract_id, ring=-1, dec_bounds=dec_bounds, ra_bounds=ra_bounds, quad=quad, is_pole=None
        )

    def summarize(self, allow_malformed: bool = False):
        """Print a summary of the converted skymap contents.

        Parameters
        ----------
        allow_malformed : bool, optional
            If True, allow malformed skymaps to be summarized.

        Raises
        ------
        ValueError
            If the skymap is malformed and allow_malformed is False.
        """
        # Check if the skymap is malformed.
        # To be well-formed, skymap must exist, and have designated file_path, metadata, tracts, and patches.
        if not allow_malformed:
            if not self.file_path.exists():
                raise ValueError(f"Skymap file does not exist: {self.file_path}")
            if not self.metadata:
                raise ValueError(f"Skymap metadata is missing: {self.file_path}")
            if not hasattr(self, "tracts"):
                raise ValueError(f"Skymap tracts are missing: {self.file_path}")
            if not hasattr(self, "patches"):
                raise ValueError(f"Skymap patches are missing: {self.file_path}")
        print("Skymap Summary")
        print("-" * 40)
        print(f"Path:               {self.file_path if self.file_path else '[unknown]'}")
        print(f"Name:               {self.metadata.get('name', '[unknown]')}")
        print(f"Generated:          {self.metadata.get('generated', '[unknown]')}")
        print(
            f"Metadata keys:      {list(self.metadata.keys()) if hasattr(self, 'metadata') else '[unknown]'}"
        )
        print(f"Number of tracts:   {self.n_tracts if hasattr(self, 'n_tracts') else '[unknown]'}")
        print(
            f"Patches per tract:  "
            f"{self.n_patches_per_tract if hasattr(self, 'n_patches_per_tract') else '[unknown]'}"
        )

        # Print file sizes for tracts and patches if available.
        if hasattr(self, "tracts_path") and self.tracts_path.exists():
            print(f"Tracts file path:   {self.tracts_path}")
            print(f"Tracts file size:   {self.tracts_path.stat().st_size / 1e6:.2f} MB")
        if hasattr(self, "patches_path") and self.patches_path.exists():
            print(f"Patches file path:  {self.patches_path}")
            print(f"Patches file size:  {self.patches_path.stat().st_size / 1e6:.2f} MB")

    def plot_patches(self, tract_patch_ids, margin=0.01, tract_outer_boundaries=None, plot_title=None):
        """Plot multiple patches in a single figure.

        Parameters
        ----------
        tract_patch_ids : list of tuples
            List of (tract_id, patch_id) tuples to plot
        margin : float, optional
            Margin (in percent) around the plotted area (default is 0.01)
        tract_outer_boundaries : int, list of int, optional
            If provided, plots the outer boundaries of specified tracts.
        plot_title : str, optional
            Title for the plot (default is None)
        """
        return plot_patches(self, tract_patch_ids, margin, tract_outer_boundaries, plot_title)


class FullVertexReader(SkymapReader):
    """Reader for full vertex format skymaps."""

    def __init__(self, file_path: str | Path):
        """Initialize the reader and load the YAML data.

        Parameters
        ----------
        file_path : str or Path
            Path to the YAML file created by FullVertexWriter
        """
        super().__init__(file_path)

        with open(self.file_path, "r") as f:
            self.data = yaml.safe_load(f)

        if "tracts" not in self.data:
            raise ValueError(f"Invalid full vertex skymap file: {file_path}")

    def get_tract(self, tract_id: int) -> TractData | None:
        """Get tract data for a given tract ID.

        Parameters
        ----------
        tract_id : int
            The tract ID to retrieve

        Returns
        -------
        TractData or None
            The tract data, or None if tract is degenerate or not found
        """
        if tract_id not in self.data["tracts"]:
            raise ValueError(f"Tract {tract_id} not found in skymap")

        content = self.data["tracts"][tract_id]
        ra_dec_vertices = content["polygon"]

        # Check for degeneracy
        unit_vecs = [UnitVector3d(LonLat.fromDegrees(ra % 360.0, dec)) for ra, dec in ra_dec_vertices]
        unique_vecs = {tuple(round(coord, 12) for coord in vec) for vec in unit_vecs}

        if len(unique_vecs) < 3:
            print(f"⚠️ Tract {tract_id} is degenerate")
            return None

        # Calculate bounds from vertices
        ras = [ra for ra, dec in ra_dec_vertices]
        decs = [dec for ra, dec in ra_dec_vertices]

        dec_bounds = (min(decs), max(decs))

        # Handle RA wrapping for bounds calculation
        ra_min, ra_max = min(ras), max(ras)
        if ra_max - ra_min > 180:  # Likely wraps around 0°
            # Find the largest gap to determine bounds
            sorted_ras = sorted(ras)
            max_gap = 0
            gap_start = 0
            for i in range(len(sorted_ras)):
                gap = (sorted_ras[(i + 1) % len(sorted_ras)] - sorted_ras[i]) % 360
                if gap > max_gap:
                    max_gap = gap
                    gap_start = sorted_ras[(i + 1) % len(sorted_ras)]
            ra_bounds = (gap_start, (gap_start - max_gap) % 360)
        else:
            ra_bounds = (ra_min, ra_max)

        return TractData(
            tract_id=tract_id,
            ring=-1,  # Full vertex format doesn't have ring info
            dec_bounds=dec_bounds,
            ra_bounds=ra_bounds,
            quad=ra_dec_vertices,  # Use the actual polygon vertices
        )

    def get_all_tracts(self) -> dict:
        """Get all tract data.

        Returns
        -------
        dict
            Dictionary mapping tract ID to TractData (or None for degenerate tracts)
        """
        tracts = {}
        for tract_id_str in self.data["tracts"]:
            tract_id = int(tract_id_str)
            tracts[tract_id] = self.get_tract(tract_id)
        return tracts

    def get_convex_polygons(self) -> dict:
        """Get ConvexPolygon objects for all tracts (legacy compatibility).

        Returns
        -------
        dict
            Dictionary mapping tract ID to ConvexPolygon or None if degenerate
        """
        poly_dict = {}

        for tract_id_str, content in self.data["tracts"].items():
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

        print(f"✅ Loaded {len(poly_dict)} tract polygons from {self.file_path}")
        return poly_dict


class RingOptimizedReader(SkymapReader):
    """A reader for ring-optimized skymaps written in YAML format.

    This class reads the YAML file and provides access to the metadata, rings, tracts, and poles.
    """

    def __init__(self, file_path: str | Path):
        """Initialize the reader and load ring-optimized skymap data.

        Parameters
        ----------
        file_path : str or Path
            Path to the YAML file created by RingOptimizedWriter
        """
        super().__init__(file_path)

        with open(self.file_path, "r") as f:
            self.data = yaml.safe_load(f)

        if "rings" not in self.data or "metadata" not in self.data:
            raise ValueError(f"Invalid ring-optimized skymap file: {file_path}")

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

    def _construct_quad(self, dec_min, dec_max, ra_start, ra_end):
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

        Returns
        -------
        list of list of float
            A list of [RA, Dec] representing the vertices of the quadrilateral.
        """
        return [
            [ra_start, dec_min],  # Lower left
            [ra_end, dec_min],  # Lower right
            [ra_end, dec_max],  # Upper right
            [ra_start, dec_max],  # Upper left
        ]

    def _construct_tract_data(self, tract_id, ring=None):
        """Construct the TractData object.

        Parameters
        ----------
        tract_id : int
            The ID of the tract.
        ring :  dict
            The ring data for the tract, if tract is not a pole.

        Returns
        -------
        TractData
            A TractData object containing the tract data.
        """
        # Out of bounds check.
        if tract_id < 0 or tract_id >= self.total_tracts:
            raise ValueError(f"Tract ID {tract_id} is out of bounds for this skymap.")

        # Handle the poles. (todo) would be nice to have this as a separate method.
        elif tract_id == 0:
            if self.poles:
                dec_min, dec_max = self.poles[0]["dec_bounds"]
                ra_start, ra_end = (
                    self.poles[0]["ra_bounds"][0],
                    self.poles[0]["ra_bounds"][1],
                )
                ra_start = ra_start % 360.0
                ra_end = ra_end % 360.0
                quad = self._construct_quad(dec_min, dec_max, ra_start, ra_end)
                return TractData(
                    tract_id=tract_id,
                    ring=-1,
                    dec_bounds=(dec_min, dec_max),
                    ra_bounds=(ra_start, ra_end),
                    quad=quad,
                )
            else:
                raise ValueError("No south pole defined in the skymap.")
        elif tract_id == self.total_tracts - 1:
            if self.poles and len(self.poles) > 1:
                dec_min, dec_max = self.poles[1]["dec_bounds"]
                ra_start, ra_end = (
                    self.poles[1]["ra_bounds"][0],
                    self.poles[1]["ra_bounds"][1],
                )
                ra_start = ra_start % 360.0
                ra_end = ra_end % 360.0
                quad = self._construct_quad(dec_min, dec_max, ra_start, ra_end)
                return TractData(
                    tract_id=tract_id,
                    ring=len(self.rings),
                    dec_bounds=(dec_min, dec_max),
                    ra_bounds=(ra_start, ra_end),
                    quad=quad,
                )
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
        return TractData(
            tract_id=tract_id,
            ring=ring["ring"],  # Extract ring number from ring dict
            dec_bounds=(dec_min, dec_max),
            ra_bounds=(ra_start, ra_end),
            quad=quad,
        )

    def get_tract(self, tract_id: int) -> TractData:
        """Get the tract data for the skymap.

        Parameters
        ----------
        tract_id : int
            The ID of the tract to retrieve.

        Returns
        -------
        TractData
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
        for ring in self.rings:
            tracts_in_ring = ring["num_tracts"]
            if tract_id <= tracts_in_ring:  # tract is in this ring
                return self._construct_tract_data(tract_id, ring)
            else:
                tract_id -= tracts_in_ring

        raise ValueError(f"Tract ID {tract_id} not found in the skymap.")

    def get_all_tracts(self) -> dict:
        """Get all tract data.

        Returns
        -------
        dict
            Dictionary mapping tract ID to TractData
        """
        tracts = {}
        for tract_id in range(self.total_tracts):
            tracts[tract_id] = self.get_tract(tract_id)
        return tracts

    def get_tract_dict(self, tract_id: int) -> dict:
        """Get tract data as dictionary (legacy compatibility).

        Parameters
        ----------
        tract_id : int
            The tract ID to retrieve

        Returns
        -------
        dict
            Dictionary representation of the tract data
        """
        tract_data = self.get_tract(tract_id)
        return tract_data.to_dict()
