import tempfile
from pathlib import Path

import numpy as np
import pytest
from skymap_convert import ConvertedSkymapReader, ConvertedSkymapWriter
from skymap_convert.test_utils import get_poly_from_tract_id, polys_are_equiv


def test_converted_skymap_structure_and_summary(lsst_skymap, capsys):
    """Test that ConvertedSkymapWriter produces readable output with correct structure."""
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir)

        writer = ConvertedSkymapWriter()
        writer.write(lsst_skymap, output_path, skymap_name="test_skymap")

        # Ensure expected files exist
        assert (output_path / "metadata.yaml").exists()
        assert (output_path / "tracts.npy").exists()
        assert (output_path / "patches.npy").exists()

        # Reader loads without error
        reader = ConvertedSkymapReader(output_path)
        assert reader.n_tracts == len(lsst_skymap)
        assert reader.n_patches_per_tract == 100
        assert reader.metadata["name"] == "test_skymap"

        # Capture and check summary output
        reader.summarize()
        output = capsys.readouterr().out
        assert "Skymap Summary" in output
        assert "test_skymap" in output


def get_patch_poly_from_ids(skymap, tract_id, patch_id):
    """Get ConvexPolygon for a given patch."""
    patch_info = skymap.generateTract(tract_id).getPatchInfo(patch_id)
    return patch_info.inner_sky_polygon


def test_converted_skymap_equivalent_to_original(lsst_skymap):
    """Test that ConvertedSkymapReader matches the original LSST SkyMap."""
    pytest.importorskip("lsst.skymap")

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir)

        writer = ConvertedSkymapWriter()
        writer.write(lsst_skymap, output_path, skymap_name="test_skymap")

        reader = ConvertedSkymapReader(output_path)

        for tract_id in range(lsst_skymap._numTracts):
            if tract_id in (0, lsst_skymap._numTracts - 1):
                continue  # Skip poles for now

            # === TRACT CHECK ===
            truth_poly = get_poly_from_tract_id(lsst_skymap, tract_id, inner=False)
            loaded_quad = reader.get_tract_vertices(tract_id)
            loaded_poly = get_poly_from_tract_id(lsst_skymap, tract_id, inner=False)  # reuse to convert quad
            loaded_poly.vertices = [np.array(v) for v in loaded_quad]

            assert polys_are_equiv(truth_poly, loaded_poly), f"Tract {tract_id} polygons not equivalent"

            # === PATCH CHECK ===
            for patch_id in range(100):
                truth_patch_poly = get_patch_poly_from_ids(lsst_skymap, tract_id, patch_id)
                loaded_patch_verts = reader.get_patch_vertices(tract_id, patch_id)

                # Convert to ConvexPolygon for comparison
                from lsst.sphgeom import ConvexPolygon, LonLat, UnitVector3d

                unit_vecs = [
                    UnitVector3d(LonLat.fromDegrees(ra % 360.0, dec)) for ra, dec in loaded_patch_verts
                ]
                loaded_patch_poly = ConvexPolygon(unit_vecs)

                assert polys_are_equiv(
                    truth_patch_poly, loaded_patch_poly
                ), f"Patch {patch_id} in tract {tract_id} not equivalent"
