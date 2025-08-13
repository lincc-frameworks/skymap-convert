from pathlib import Path

import pytest
from skymap_convert import ConvertedSkymapReader
from skymap_convert.utils import load_pickle_skymap

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
TEST_DIR = PACKAGE_ROOT / "tests"
RAW_SKYMAP_DIR = TEST_DIR / "data" / "raw_skymaps"


def pytest_addoption(parser):
    """Add command line options for pytest.

    Parameters
    ----------
    parser : pytest.Parser
        The pytest argument parser to add options to.
    """
    parser.addoption(
        "--longrun", action="store_true", dest="longrun", default=False, help="enable longrun decorated tests"
    )


def pytest_configure(config):
    """Configure pytest to skip longrun tests unless --longrun is specified.

    Parameters
    ----------
    config : pytest.Config
        The pytest configuration object.

    Notes
    -----
    This function automatically excludes tests marked with @pytest.mark.longrun
    unless the --longrun flag is explicitly provided on the command line.
    """
    # If the --longrun option is not specified, skip tests marked with @pytest.mark.longrun
    if not config.option.longrun:
        config.option.markexpr = "not longrun"


@pytest.fixture(scope="session")
def converted_skymap_reader():
    """Fixture that provides a ConvertedSkymapReader for testing.

    Returns
    -------
    ConvertedSkymapReader
        A reader instance using the "lsst_skymap" preset.

    Note
    ----
    This uses the pre-converted LSST skymap for performance reasons.
    """

    return ConvertedSkymapReader(preset="lsst_skymap")


@pytest.fixture(scope="session")
def lsst_skymap():
    """Fixture to provide an LSST skymap object for testing.

    Returns
    -------
    lsst.skymap.SkyMap
        The loaded LSST skymap object from the test data

    Notes
    -----
    This fixture automatically skips tests if the lsst.skymap package
    is not available. The skymap is loaded from a pickle file in the
    test data directory.
    """
    pytest.importorskip("lsst.skymap")
    return load_pickle_skymap(RAW_SKYMAP_DIR / "skyMap_lsst_cells_v1_skymaps.pickle")
