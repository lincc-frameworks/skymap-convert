import pytest
import yaml
from skymap_convert.utils import (
    get_quad_from_patch_id,
    get_quad_from_tract_id,
    quads_are_equiv,
)
from tqdm import tqdm

TRACT_SAMPLES = [1, 250, 1899]
PATCH_SAMPLES = [0, 42, 99]


@pytest.mark.longrun
def test_writer_creates_expected_files(lsst_skymap, written_skymap_data):
    """Test that ConvertedSkymapWriter creates all expected output files.

    Parameters
    ----------
    lsst_skymap : lsst.skymap.SkyMap
        The original LSST skymap object to write.
    written_skymap_data : dict
        Fixture providing written skymap data with output_dir, skymap_name, and reader.

    Notes
    -----
    Verifies that:
    - All expected files (metadata.yaml, tracts.npy, patches.npy.gz) are created
    - No temporary files are left behind
    - Files have non-zero size
    """
    pytest.importorskip("lsst.skymap")
    pytest.importorskip("lsst.sphgeom")

    output_dir = written_skymap_data["output_dir"]
    skymap_name = written_skymap_data["skymap_name"]

    # Check expected files exist
    assert (output_dir / "metadata.yaml").exists()
    assert (output_dir / "tracts.npy").exists()
    assert (output_dir / "patches.npy.gz").exists()

    # Check temporary files are cleaned up
    assert not (output_dir / "patches.npy").exists()

    # Check files have content
    assert (output_dir / "metadata.yaml").stat().st_size > 0
    assert (output_dir / "tracts.npy").stat().st_size > 0
    assert (output_dir / "patches.npy.gz").stat().st_size > 0

    # Load and check metadata
    with open(output_dir / "metadata.yaml", "r") as f:
        metadata = yaml.safe_load(f)

    assert metadata["name"] == skymap_name
    assert metadata["n_tracts"] == len(lsst_skymap)
    assert metadata["n_patches_per_tract"] == 100
    assert metadata["format_version"] == 1
    assert "generated" in metadata
    assert metadata["generated"].endswith("Z")  # ISO format with UTC


@pytest.mark.longrun
def test_writer_sample_tracts(lsst_skymap, written_skymap_data):
    """Test complete round-trip: write with ConvertedSkymapWriter, read with ConvertedSkymapReader.

    Parameters
    ----------
    lsst_skymap : lsst.skymap.SkyMap
        The original LSST skymap object.
    written_skymap_data : dict
        Fixture providing written skymap data with output_dir, skymap_name, and reader.

    Notes
    -----
    Verifies that data written by ConvertedSkymapWriter can be correctly
    read by ConvertedSkymapReader and produces identical results.
    """
    pytest.importorskip("lsst.skymap")
    pytest.importorskip("lsst.sphgeom")

    reader = written_skymap_data["reader"]
    skymap_name = written_skymap_data["skymap_name"]

    # Verify metadata matches
    assert reader.metadata["name"] == skymap_name
    assert reader.n_tracts == len(lsst_skymap)
    assert reader.n_patches_per_tract == 100

    # Verify a few sample geometries match original
    for tract_id in TRACT_SAMPLES:
        original_tract = get_quad_from_tract_id(lsst_skymap, tract_id, inner=True)
        round_trip_tract = reader.get_tract_vertices(tract_id)

        assert quads_are_equiv(original_tract, round_trip_tract), f"Round-trip failed for tract {tract_id}"

        for patch_id in PATCH_SAMPLES:
            original_patch = get_quad_from_patch_id(lsst_skymap, tract_id, patch_id)
            round_trip_patch = reader.get_patch_vertices(tract_id, patch_id)

            assert quads_are_equiv(
                original_patch, round_trip_patch
            ), f"Round-trip failed for patch {patch_id} in tract {tract_id}"


@pytest.mark.longrun
def test_writer_full_skymap_integrity(lsst_skymap, written_skymap_data):
    """Comprehensive test that writer preserves all tract and patch geometries.

    Parameters
    ----------
    lsst_skymap : lsst.skymap.SkyMap
        The original LSST skymap object.
    written_skymap_data : dict
        Fixture providing written skymap data with output_dir, skymap_name, and reader.

    Notes
    -----
    This is a comprehensive test that verifies every tract and patch geometry
    is preserved correctly using the session-scoped written skymap data.
    """
    if True:
        pytest.skip(
            "Skipping full integrity test as it takes an exceptionally long time to run. "
            "Recommend running periodically, especially after major changes to the writer. "
            "To run, manually remove the branching logic in the code of the test."
        )
    else:
        pytest.importorskip("lsst.skymap")
        pytest.importorskip("lsst.sphgeom")

        reader = written_skymap_data["reader"]

        tract_ids = range(len(lsst_skymap))
        for tract_id in tqdm(tract_ids, desc="Verifying written tracts", leave=False):
            # Verify tract geometry
            truth_quad = get_quad_from_tract_id(lsst_skymap, tract_id, inner=True)
            written_quad = reader.get_tract_vertices(tract_id)

            assert quads_are_equiv(truth_quad, written_quad), f"Tract {tract_id} geometry not preserved"

            # Verify all patch geometries
            for patch_id in range(100):
                truth_patch = get_quad_from_patch_id(lsst_skymap, tract_id, patch_id)
                written_patch = reader.get_patch_vertices(tract_id, patch_id)

                assert quads_are_equiv(
                    truth_patch, written_patch
                ), f"Patch {patch_id} in tract {tract_id} geometry not preserved"
