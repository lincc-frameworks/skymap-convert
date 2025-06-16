import pytest
from skymap_convert import load_polygons_ra_dec, write_polygons_ra_dec
from skymap_convert.test_utils import get_poly_from_tract_id, polys_are_equiv
import tempfile
from pathlib import Path


@pytest.mark.parametrize("inner", [True, False])
def test_loaded_polygons_equivalent(lsst_skymap, inner):
    """Test that polygons written and reloaded from disk are equivalent to ground truth."""

    # use tmpdir fixture to create a temporary directory for the test
    with tempfile.TemporaryDirectory() as tmpdir:
        converted_skymap_path = Path(tmpdir) / ("inner_polygons.yaml" if inner else "outer_polygons.yaml")

        # Write the skymap to disk
        write_polygons_ra_dec(lsst_skymap, converted_skymap_path, inner=inner)

        # Load from disk
        loaded_poly_map = load_polygons_ra_dec(converted_skymap_path)

    for tract_id in range(lsst_skymap._numTracts):
        ground_truth = get_poly_from_tract_id(lsst_skymap, tract_id, inner=inner)
        loaded = loaded_poly_map.get(tract_id)

        if tract_id == 0 or tract_id == lsst_skymap._numTracts - 1:
            # Skip first and last tract for now, but todo
            continue
        else:
            assert loaded is not None, f"Tract {tract_id} missing from loaded polygons."
            assert polys_are_equiv(ground_truth, loaded), f"Tract {tract_id} polygons are not equivalent"
