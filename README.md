<!--
[![PyPI](https://img.shields.io/pypi/v/skymap-convert?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/skymap_convert/)
[![Codecov](https://codecov.io/gh/lincc-frameworks/skymap-convert/branch/main/graph/badge.svg)](https://codecov.io/gh/lincc-frameworks/skymap-convert)
-->
[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/lincc-frameworks/skymap-convert/smoke-test.yml)](https://github.com/lincc-frameworks/skymap-convert/actions/workflows/smoke-test.yml)

# Skymap Convert
A dependency-light package for skymap operations--use the LSST skymap without the LSST stack.

## Quick start

1. Install the package.
```bash
git clone https://github.com/lincc-frameworks/skymap-convert.git
cd skymap-convert
pip install .
```

2. Import a skymap reader.
```python
import skymap_convert

reader = skymap_convert.ConvertedSkymapReader(preset="lsst_skymap")
```

3. Optionally, call `summarize` to take a peek at the convents of the converted skymap.
```python
reader.summarize()

# Prints:
# 
# Skymap Summary
# ----------------------------------------
# Path:               /path/to/skymap-convert/converted_skymaps/lsst_skymap
# Name:               lsst_skymap
# Generated:          2025-07-01T18:11:20.873149Z
# ...
```
4. Use `get_tract_vertices` and `get_patch_vertices` to access the data.  
*(See [Reading and using a converted skymap.ipynb](https://github.com/lincc-frameworks/skymap-convert/blob/main/docs/notebooks/02%20-%20Reading%20and%20using%20a%20converted%20skymap.ipynb) for more details.)*
```python
reader.get_tract_vertices(42)

# Returns:
#
# [[86.16241626810654, -88.63764259611838],
# [92.73276185933494, -88.5876043887882],
# [90.57872844947276, -88.43062126829582],
# [84.63000467541433, -88.47549635501055]]
```
```python
reader.get_patch_vertices(1234, 56)

# Returns:
#
# [[353.19436746094334, -61.733705740129906],
# [353.5462936615678, -61.73505840062234],
# [353.54818789227375, -61.568395336612376],
# [353.19815518019567, -61.567052069980534]]
```

5. To plot, call `plot_patches`.  
*(See [Skymap plotting methods.ipynb](https://github.com/lincc-frameworks/skymap-convert/blob/main/docs/notebooks/Skymap%20plotting%20methods.ipynb) for more details.)*
```python
reader.plot_patches(
    [
        (60, 0),
        (61, 8)
    ],
    tract_outer_boundaries=60
)
```

## For more information
### A general overview on skymaps and how they are structured
See [Skymaps overview](https://github.com/lincc-frameworks/skymap-convert/blob/main/docs/notebooks/00%20-%20Skymaps%20overview.ipynb).

### Converting your own skymaps
 See [Converting from LSST skymaps.ipynb](https://github.com/lincc-frameworks/skymap-convert/blob/main/docs/notebooks/01%20-%20Converting%20from%20LSST%20skymaps.ipynb).


