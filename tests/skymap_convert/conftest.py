from pathlib import Path

import pytest
from skymap_convert import ConvertedSkymapReader, ConvertedSkymapWriter
from skymap_convert.utils import load_pickle_skymap

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
TEST_DIR = PACKAGE_ROOT / "tests"
RAW_SKYMAP_DIR = TEST_DIR / "data" / "raw_skymaps"
SKYMAP_OUT_DIR = PACKAGE_ROOT / "converted_skymaps"


def pytest_addoption(parser):
    """Add command line options for pytest."""
    parser.addoption(
        "--longrun", action="store_true", dest="longrun", default=False, help="enable longrun decorated tests"
    )


def pytest_configure(config):
    """Configure pytest to skip longrun tests unless --longrun is specified."""
    # If the --longrun option is not specified, skip tests marked with @pytest.mark.longrun
    if not config.option.longrun:
        config.option.markexpr = "not longrun"


@pytest.fixture(scope="session")
def converted_skymap_reader(tmp_path_factory, lsst_skymap):
    """Fixture that writes the converted skymap once and returns a reader."""
    tmp_dir = tmp_path_factory.mktemp("converted_skymap")
    output_path = tmp_dir / "converted"

    print("Writing converted skymap for the testing session...")

    writer = ConvertedSkymapWriter()
    writer.write(lsst_skymap, output_path, skymap_name="test_skymap")

    print(f"Converted skymap written to {output_path}")

    return ConvertedSkymapReader(output_path)


@pytest.fixture(scope="session")
def skymap_out_dir():
    """Fixture to provide the output directory for skymap polygons."""
    return SKYMAP_OUT_DIR


@pytest.fixture(scope="session")
def lsst_skymap():
    """Fixture to provide a LSST skymap object."""
    pytest.importorskip("lsst.skymap")
    return load_pickle_skymap(RAW_SKYMAP_DIR / "skyMap_lsst_cells_v1_skymaps.pickle")
