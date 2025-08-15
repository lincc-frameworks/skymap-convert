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
    The reader will be automatically cleaned up at the end of the test session.
    """
    reader = ConvertedSkymapReader(preset="lsst_skymap")

    # Yield the reader for use in tests
    yield reader

    # Cleanup after all tests are done
    reader.cleanup()


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


@pytest.fixture(scope="session")
def written_skymap_data(lsst_skymap, tmp_path_factory, request):
    """Session-scoped fixture that writes skymap data once and provides paths for testing.

    Parameters
    ----------
    lsst_skymap : lsst.skymap.SkyMap
        The original LSST skymap object to write.
    tmp_path_factory : pytest.TempPathFactory
        Pytest factory for creating temporary directories.
    request : pytest.FixtureRequest
        Pytest request object to access command line options.

    Returns
    -------
    dict
        Dictionary containing:
        - 'output_dir': Path to the directory containing written skymap files
        - 'skymap_name': Name used for the skymap
        - 'reader': ConvertedSkymapReader instance for the written data

    Notes
    -----
    This fixture only runs if --longrun is specified, otherwise it skips.
    The expensive skymap writing operation happens only once per test session.
    The reader will be automatically cleaned up at the end of the test session.
    """
    # Skip if not running longrun tests
    if not request.config.option.longrun:
        pytest.skip("Skipping written_skymap_data fixture - requires --longrun")

    pytest.importorskip("lsst.skymap")
    pytest.importorskip("lsst.sphgeom")

    from skymap_convert.skymap_writers import ConvertedSkymapWriter

    # Create session-scoped temporary directory
    output_dir = tmp_path_factory.mktemp("session_skymap_data")
    skymap_name = "test_session_skymap"

    # Write the skymap once
    writer = ConvertedSkymapWriter()
    writer.write(lsst_skymap, output_dir, skymap_name)

    # Create reader for the written data
    reader = ConvertedSkymapReader(output_dir)

    data = {"output_dir": output_dir, "skymap_name": skymap_name, "reader": reader}

    # Yield the data for use in tests
    yield data

    # Cleanup after all tests are done
    reader.cleanup()
