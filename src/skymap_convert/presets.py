"""
Built-in skymap presets and paths.

This module provides easy access to the converted skymaps that come bundled
with the skymap-convert package, eliminating the need to manually construct
paths to the built-in skymap data.
"""

# import sys
from pathlib import Path

# Get the package root directory
_PACKAGE_ROOT = Path(__file__).parent.parent.parent
_CONVERTED_SKYMAPS_DIR = _PACKAGE_ROOT / "converted_skymaps"


def get_preset_path(preset_name: str) -> Path:
    """Get the path to a built-in skymap preset.

    Parameters
    ----------
    preset_name : str
        Name of the preset skymap

    Returns
    -------
    Path
        Path to the preset skymap directory

    Raises
    ------
    ValueError
        If the preset name is not recognized
    FileNotFoundError
        If the preset directory doesn't exist
    """
    preset_path = _CONVERTED_SKYMAPS_DIR / preset_name

    if not preset_path.exists():
        available = list_available_presets()
        raise FileNotFoundError(f"Preset '{preset_name}' not found. Available presets: {available}")

    return preset_path


def list_available_presets() -> list[str]:
    """List all available built-in skymap presets.

    Returns
    -------
    list[str]
        List of available preset names
    """
    if not _CONVERTED_SKYMAPS_DIR.exists():
        return []

    presets = []
    for item in _CONVERTED_SKYMAPS_DIR.iterdir():
        if item.is_dir() and not item.name.startswith("."):
            presets.append(item.name)

    return sorted(presets)


def get_preset_info() -> dict[str, dict[str, str]]:
    """Get information about available presets.

    Returns
    -------
    dict[str, dict[str, str]]
        Dictionary mapping preset names to their metadata
    """
    info = {}

    for preset_name in list_available_presets():
        preset_path = _CONVERTED_SKYMAPS_DIR / preset_name
        metadata_path = preset_path / "metadata.yaml"

        if metadata_path.exists():
            import yaml

            with open(metadata_path, "r") as f:
                metadata = yaml.safe_load(f)
                info[preset_name] = {
                    "path": str(preset_path),
                    "name": metadata.get("name", "Unknown"),
                    "generated": metadata.get("generated", "Unknown"),
                    "n_tracts": metadata.get("n_tracts", "Unknown"),
                    "n_patches_per_tract": metadata.get("n_patches_per_tract", "Unknown"),
                }
        else:
            info[preset_name] = {
                "path": str(preset_path),
                "name": preset_name,
                "generated": "Unknown",
                "n_tracts": "Unknown",
                "n_patches_per_tract": "Unknown",
            }

    return info


# class _PresetsModule:
#     """Helper class to provide attribute access to preset paths."""

#     @property
#     def lsst_skymap(self) -> Path:
#         """Path to the built-in LSST skymap preset."""
#         return get_preset_path("lsst_skymap")

#     def get_preset_path(self, preset_name: str) -> Path:
#         """Get the path to a built-in skymap preset."""
#         return get_preset_path(preset_name)

#     def list_available_presets(self) -> list[str]:
#         """List all available built-in skymap presets."""
#         return list_available_presets()

#     def get_preset_info(self) -> dict[str, dict[str, str]]:
#         """Get information about available presets."""
#         return get_preset_info()


# # Replace the module with our custom class instance
# sys.modules[__name__] = _PresetsModule()
