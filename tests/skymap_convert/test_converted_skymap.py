import pytest
from skymap_convert.test_utils import (
    get_quad_from_patch_id,
    get_quad_from_tract_id,
    quads_are_equiv,
)
from tqdm import tqdm

TRACT_SAMPLES = [1, 250, 1899]
PATCH_SAMPLES = [0, 42, 99]


def test_converted_skymap_structure_and_summary(lsst_skymap, converted_skymap_reader, capsys):
    """Test that ConvertedSkymapWriter produces readable output with correct structure."""
    # Ensure expected files exist
    assert (converted_skymap_reader.metadata_path).exists()
    assert (converted_skymap_reader.tracts_path).exists()
    assert (converted_skymap_reader.patches_path).exists()

    # Reader loads without error
    assert converted_skymap_reader.n_tracts == len(lsst_skymap)
    assert converted_skymap_reader.n_patches_per_tract == 100
    assert converted_skymap_reader.metadata["name"] == "test_skymap"

    # Capture and check summary output
    converted_skymap_reader.summarize()
    output = capsys.readouterr().out
    assert "Skymap Summary" in output
    assert "test_skymap" in output


@pytest.mark.parametrize("tract_id", TRACT_SAMPLES)
def test_sample_tracts_and_patches(lsst_skymap, converted_skymap_reader, tract_id, tmp_path):
    """Quick check that a subset of tracts and patches match after conversion."""
    pytest.importorskip("lsst.skymap")

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
    """Test that ConvertedSkymapReader matches the original LSST SkyMap."""
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
