from pathlib import Path

import pytest

from skymap_convert.io import load_pickle_skymap

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
TEST_DIR = PACKAGE_ROOT / "tests"
RAW_SKYMAP_DIR = TEST_DIR / "data" / "raw_skymaps"
SKYMAP_OUT_DIR = PACKAGE_ROOT / "converted_skymaps"


@pytest.fixture
def skymap_out_dir():
    """Fixture to provide the output directory for skymap polygons."""
    return SKYMAP_OUT_DIR


@pytest.fixture
def lsst_skymap():
    """Fixture to provide a LSST skymap object."""
    return load_pickle_skymap(RAW_SKYMAP_DIR / "skyMap_lsst_cells_v1_skymaps.pickle")
