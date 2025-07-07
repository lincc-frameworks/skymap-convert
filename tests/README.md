# Notes on included tests

## include or skip
A lot of our tests use `pytest.importorskip("lsst.skymap")`. These will be skipped in environments lacking the `lsst.skymap` package, ie, environments that do not have access to the LSST stack. Writing a new converted skymap will always require access to the LSST stack, but reading a converted skymap in general should not (unless checking against a `lsst.skymap` type), so this allows users to still run any reading-related tests.

## longrun test(s) and tqdm
Use `pytest -k test_converted_skymap_equivalent_to_original -s --longrun` to run the test that fully checks the converted skymap (each tract, each patch in each tract, and each vertex in each patch) against the original lsst.skymap type skymap. It takes a little while to run, so the `logrun` decorator has been added, and this test will automatically be skipped on a simple call to `pytest`.

Note that the '-s' flag is necessary for tqdm to show progress through each tract for the latter part of that test run.

## session-scoped converted skymap reader
Also, note that we do execute the writer at the beginning of a testing session (unless you're specifically testing a single test that does not use the `converted_skymap_reader` reader), instead of checking for any existing converted skymap files. However, this also takes a while, so we have a `scope="session"` decorator to only do this once per test session. 