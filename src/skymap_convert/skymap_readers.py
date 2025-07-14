import gzip
import shutil
import tempfile
from pathlib import Path

import numpy as np
import yaml

from .plotting import plot_patches


class ConvertedSkymapReader:
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
        # Use preset if provided, otherwise check for file_path
        if preset is not None:
            from .presets import get_preset_path

            file_path = get_preset_path(preset)
        elif file_path is None:
            raise ValueError("Either file_path or preset must be specified")

        # Set the file path
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"Skymap file not found: {file_path}")
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
