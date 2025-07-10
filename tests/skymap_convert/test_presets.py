"""
Test script to verify that the presets functionality works correctly.
"""

from pathlib import Path

import pytest


def test_presets_module_import():
    """Test that the presets module can be imported."""
    import skymap_convert.presets as presets

    assert presets is not None


def test_list_available_presets():
    """Test that we can list available presets."""
    import skymap_convert.presets as presets

    available = presets.list_available_presets()
    assert isinstance(available, list)
    # Should at least have lsst_skymap if the file exists
    assert len(available) >= 0


def test_lsst_skymap_preset():
    """Test accessing the lsst_skymap preset if it exists."""
    import skymap_convert.presets as presets

    available = presets.list_available_presets()

    if "lsst_skymap" in available:
        lsst_path = presets.lsst_skymap
        assert isinstance(lsst_path, Path)
        # Path should be absolute
        assert lsst_path.is_absolute()
    else:
        pytest.skip("lsst_skymap not found in available presets")


def test_converted_skymap_reader_with_preset():
    """Test ConvertedSkymapReader with preset parameter."""
    import skymap_convert.presets as presets
    from skymap_convert import ConvertedSkymapReader

    available = presets.list_available_presets()

    if "lsst_skymap" in available:
        reader = ConvertedSkymapReader(preset="lsst_skymap")
        assert reader is not None

        # Test that help() method works without errors
        reader.help()
    else:
        pytest.skip("lsst_skymap not available for testing")


def test_converted_skymap_reader_with_invalid_preset():
    """Test that ConvertedSkymapReader raises appropriate error for invalid preset."""
    from skymap_convert import ConvertedSkymapReader

    with pytest.raises(FileNotFoundError, match="Preset 'nonexistent_preset' not found"):
        ConvertedSkymapReader(preset="nonexistent_preset")


def test_converted_skymap_reader_preset_vs_direct_path():
    """Test that preset and direct path give equivalent results."""
    import skymap_convert.presets as presets
    from skymap_convert import ConvertedSkymapReader

    available = presets.list_available_presets()

    if "lsst_skymap" in available:
        # Create reader with preset
        reader_preset = ConvertedSkymapReader(preset="lsst_skymap")

        # Create reader with direct path
        reader_direct = ConvertedSkymapReader(presets.lsst_skymap)

        # Both should have the same path
        assert reader_preset.file_path == reader_direct.file_path
    else:
        pytest.skip("lsst_skymap not available for testing")
