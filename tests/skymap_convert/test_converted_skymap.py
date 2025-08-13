import pytest
from skymap_convert.utils import (
    get_quad_from_patch_id,
    get_quad_from_tract_id,
    quads_are_equiv,
)
from tqdm import tqdm

TRACT_SAMPLES = [1, 250, 1899]
PATCH_SAMPLES = [0, 42, 99]


def test_converted_skymap_structure_and_summary(lsst_skymap, converted_skymap_reader, capsys):
    """Test that ConvertedSkymapReader produces readable output with correct structure.

    Parameters
    ----------
    lsst_skymap : lsst.skymap.SkyMap
        The original LSST skymap object for comparison.
    converted_skymap_reader : ConvertedSkymapReader
        Reader instance for the converted skymap.
    capsys : pytest.CaptureFixture
        Pytest fixture for capturing stdout/stderr.

    Notes
    -----
    Verifies that:
    - All expected files (metadata.yaml, tracts.npy, patches.npy) exist
    - Reader loads metadata correctly
    - Summary output contains expected content
    """
    # Ensure expected files exist
    assert (converted_skymap_reader.metadata_path).exists()
    assert (converted_skymap_reader.tracts_path).exists()
    assert (converted_skymap_reader.patches_path).exists()

    # Reader loads without error
    assert converted_skymap_reader.n_tracts == len(lsst_skymap)
    assert converted_skymap_reader.n_patches_per_tract == 100
    assert converted_skymap_reader.metadata["name"] == "lsst_skymap"

    # Capture and check summary output
    converted_skymap_reader.summarize()
    output = capsys.readouterr().out
    assert "Skymap Summary" in output
    assert "lsst_skymap" in output


@pytest.mark.parametrize("tract_id", TRACT_SAMPLES)
def test_sample_tracts_and_patches(lsst_skymap, converted_skymap_reader, tract_id, tmp_path):
    """Verify that sample tracts and patches match between original and converted skymaps.

    Parameters
    ----------
    lsst_skymap : lsst.skymap.SkyMap
        The original LSST skymap object for comparison.
    converted_skymap_reader : ConvertedSkymapReader
        Reader instance for the converted skymap.
    tract_id : int
        The tract ID to test (parametrized from TRACT_SAMPLES).
    tmp_path : Path
        Pytest fixture providing temporary directory path.

    Notes
    -----
    Tests a representative subset of tracts and patches for performance.
    For each sample tract, compares:
    - Tract boundary vertices
    - Sample patch vertices (patches 0, 42, 99)

    Requires lsst.skymap and lsst.sphgeom packages for reading and processing the original skymap.
    If these packages are not available, the test will be skipped.
    """
    pytest.importorskip("lsst.skymap")
    pytest.importorskip("lsst.sphgeom")

    # Tract comparison
    truth_quad = get_quad_from_tract_id(lsst_skymap, tract_id, inner=True)
    loaded_quad = converted_skymap_reader.get_tract_vertices(tract_id)

    assert quads_are_equiv(truth_quad, loaded_quad), f"Tract {tract_id} mismatch"

    # Patch comparisons
    for patch_id in PATCH_SAMPLES:
        truth_patch = get_quad_from_patch_id(lsst_skymap, tract_id, patch_id)
        loaded_patch = converted_skymap_reader.get_patch_vertices(tract_id, patch_id)

        assert quads_are_equiv(truth_patch, loaded_patch), f"Patch {patch_id} in tract {tract_id} mismatch"


@pytest.mark.longrun
def test_converted_skymap_equivalent_to_original(lsst_skymap, converted_skymap_reader):
    """Comprehensive test that ConvertedSkymapReader exactly matches the original LSST SkyMap.

    Parameters
    ----------
    lsst_skymap : lsst.skymap.SkyMap
        The original LSST skymap object for comparison
    converted_skymap_reader : ConvertedSkymapReader
        Reader instance for the converted skymap

    Notes
    -----
    This is a comprehensive test that compares every tract and patch between
    the original and converted skymaps. It is marked as @pytest.mark.longrun
    because it can take several minutes to complete.

    Run with: pytest --longrun to include this test.

    For each tract in the skymap:
    - Compares tract boundary vertices
    - Compares all 100 patch vertices within the tract

    Requires lsst.skymap and lsst.sphgeom packages.
    """
    pytest.importorskip("lsst.skymap")

    tract_ids = range(lsst_skymap._numTracts)
    for tract_id in tqdm(tract_ids, desc="Checking tracts", leave=False):
        # Tract check
        truth_quad = get_quad_from_tract_id(lsst_skymap, tract_id, inner=True)
        loaded_quad = converted_skymap_reader.get_tract_vertices(tract_id)

        assert quads_are_equiv(truth_quad, loaded_quad), f"Tract {tract_id} quads not equivalent"

        # Patch check
        for patch_id in range(100):
            truth_patch_quad = get_quad_from_patch_id(lsst_skymap, tract_id, patch_id)
            loaded_patch_quad = converted_skymap_reader.get_patch_vertices(tract_id, patch_id)

            assert quads_are_equiv(
                truth_patch_quad, loaded_patch_quad
            ), f"Patch {patch_id} in tract {tract_id} not equivalent"
