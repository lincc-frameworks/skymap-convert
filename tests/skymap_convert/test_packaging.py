import skymap_convert


def test_version():
    """Check to see that we can get the package version"""
    assert skymap_convert.__version__ is not None
