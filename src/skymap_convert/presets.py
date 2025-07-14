"""
Built-in skymap presets and paths.

This module provides easy access to the converted skymaps that come bundled
with the skymap-convert package, eliminating the need to manually construct
paths to the built-in skymap data.
"""

from importlib.resources import files
from pathlib import Path


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
    presets = files("skymap_convert.converted_skymaps")
    preset_path = presets / preset_name

    try:
        # Check if we can iterate the directory (this tests existence)
        list(preset_path.iterdir())
        return Path(preset_path)
    except (OSError, FileNotFoundError) as err:
        available = list_available_presets()
        raise FileNotFoundError(f"Preset '{preset_name}' not found. Available presets: {available}") from err


def list_available_presets() -> list[str]:
    """List all available built-in skymap presets.

    Returns
    -------
    list[str]
        List of available preset names
    """
    try:
        presets_dir = files("skymap_convert.converted_skymaps")

        preset_names = []
        for item in presets_dir.iterdir():
            if item.is_dir() and not item.name.startswith("."):
                preset_names.append(item.name)

        return sorted(preset_names)
    except (FileNotFoundError, ModuleNotFoundError):
        return []


def get_preset_info() -> dict[str, dict[str, str]]:
    """Get information about available presets.

    Returns
    -------
    dict[str, dict[str, str]]
        Dictionary mapping preset names to their metadata
    """
    info = {}

    for preset_name in list_available_presets():
        preset_path = get_preset_path(preset_name)
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
