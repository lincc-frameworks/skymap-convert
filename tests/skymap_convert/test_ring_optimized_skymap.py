import tempfile
from pathlib import Path

import pytest
from skymap_convert import RingOptimizedReader, RingOptimizedWriter
from skymap_convert.test_utils import get_poly_from_tract_id, polys_are_equiv


@pytest.mark.parametrize("inner", [True])
def test_ring_optimized_skymap_read_write(lsst_skymap, inner):
    """Test that polygons calculated by data written and reloaded from disk are equivalent to ground truth."""
    pytest.importorskip("lsst.skymap")

    # use tmpdir fixture to create a temporary directory for the test
    with tempfile.TemporaryDirectory() as tmpdir:
        converted_skymap_path = Path(tmpdir) / ("inner_polygons.yaml" if inner else "outer_polygons.yaml")

        # Write the skymap to disk using OOP writer
        writer = RingOptimizedWriter()
        writer.write(lsst_skymap, converted_skymap_path, inner=inner)

        # Load from disk using OOP reader
        reader = RingOptimizedReader(converted_skymap_path)

    for tract_id in range(lsst_skymap._numTracts):
        ground_truth = get_poly_from_tract_id(lsst_skymap, tract_id, inner=inner)

        if tract_id == 0 or tract_id == lsst_skymap._numTracts - 1:
            # Skip first and last tract for now, but todo
            continue

        # Get tract data and convert to ConvexPolygon
        tract_data = reader.get_tract(tract_id)
        loaded = None if tract_data is None else tract_data.to_convex_polygon()

        assert loaded is not None, f"Tract {tract_id} missing from loaded polygons."
        assert polys_are_equiv(
            ground_truth, loaded, atol=1e-8
        ), f"Tract {tract_id} polygons are not equivalent"
