import pytest
from skymap_to_poly_coords import load_polygons_ra_dec
from skymap_to_poly_coords.test_utils import polys_are_equiv, get_poly_from_tract_id


@pytest.mark.parametrize("inner", [True, False])
def test_loaded_polygons_equivalent(lsst_skymap, skymap_out_dir, inner):
    """Test that polygons written and reloaded from disk are equivalent to ground truth."""

    # Choose path based on inner/outer
    yaml_path = skymap_out_dir / ("inner_polys.yaml" if inner else "outer_polys.yaml")

    # Load from disk
    loaded_poly_map = load_polygons_ra_dec(yaml_path)

    for tract_id in range(lsst_skymap._numTracts):
        ground_truth = get_poly_from_tract_id(lsst_skymap, tract_id, inner=inner)
        loaded = loaded_poly_map.get(tract_id)

        if tract_id == 0 or tract_id == lsst_skymap._numTracts - 1:
            # Skip first and last tract for now, but todo
            continue
        else:
            assert loaded is not None, f"Tract {tract_id} missing from loaded polygons."
            assert polys_are_equiv(ground_truth, loaded), f"Tract {tract_id} polygons are not equivalent"
