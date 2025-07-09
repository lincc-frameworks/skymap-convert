# Skymap Convert

<!--
[![PyPI](https://img.shields.io/pypi/v/skymap-convert?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/skymap_convert/)
[![Codecov](https://codecov.io/gh/lincc-frameworks/skymap-convert/branch/main/graph/badge.svg)](https://codecov.io/gh/lincc-frameworks/skymap-convert)
-->
[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/lincc-frameworks/skymap-convert/smoke-test.yml)](https://github.com/lincc-frameworks/skymap-convert/actions/workflows/smoke-test.yml)


A dependency-light package for working with skymaps--use the LSST skymap without needing the LSST stack.

## Quick start


1. Install the package
```bash
pip install skymap-convert
```

2. Import a skymap reader
```python
import skymap-convert
reader = skymap_convert.ConvertedSkymapReader(converted_skymap_path)
```

3. Optionally, call `summarize` to take a peek at the convents of the converted skymap
```python
reader.summarize()
```

4. Use `get_tract_vertices` and `get_patch_vertices` to access the data
```python
reader.get_tract_vertices(42)
reader.get_patch_vertices(1234, 56)
```

5. To plot, call `plot_patches`  
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
