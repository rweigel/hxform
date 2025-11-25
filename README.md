# `hxform`

Heliophysics reference frame ("coordinate") transforms in Python using various libraries, Python packages, and web services using a common interface.

Libraries include
[`cxform`](https://github.com/edsantiago/cxform)
and
[`Geopack-08`](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/geopack/).
Packages include 
[`PySPEDAS`](https://github.com/spedas/pyspedas),
[`SpacePy`](https://github.com/spacepy/spacepy),
[`SpiceyPy`](https://github.com/AndrewAnnex/SpiceyPy),
and
[`SunPy`](https://github.com/sunpy/sunpy).
Web services include
[`SSCWeb Coordinate Calculator`](https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi).

# Motivation

Over time, we have needed to use various coordinate transform libraries in Python. Multiple libraries have been used because

1. no single library had all needed transforms,
2. some libraries were too slow for an application,
3. some libraries had errors for certain transforms, and
4. some libraries did not interface with ParaView's Python.

To facilitate comparison, we have developed a single interface. To transform a vector `v` at time `t` from frame  `frame_in` to `frame_out`, using library `LIB`, use

```python
vt = hxform.transform(v, t, frame_in, frame_out, lib=LIB)
```

See the [`demo/`](https://github.com/rweigel/hxform/tree/main/demo) directory for example usage.

# Overview

The list of supported libraries and frames is in [`info.py`](https://github.com/rweigel/hxform/blob/main/hxform/info.py).

# Install

Tested with Python 3.{11,12,13,14}.

```bash
git clone https://github.com/rweigel/hxform
cd hxform
pip install -e .
```

To test if the installation was successful, execute

```
python demo/basic.py
```

# Tests

## Single version of Python

See the files in the [test](https://github.com/rweigel/hxform/tree/master/test) directory. The result from executing the files is stored in a `.log` file.

To run, execute

```bash
pip install pytest
pytest # Executes files in ./test
```

## Multiple versions of Python

```
python -m pip install --upgrade 'nox[pbs]'
nox -s tests
# To test using only a single version of Python
nox -s tests --python 3.12
```

# Development

## Native interfaces

Examples of using the native interface for other Python packages are in the [`demo-native`](https://github.com/rweigel/hxform/tree/main/demo-native) directory.

## Geopack-08 and cxform

[cxform (C)](https://github.com/edsantiago/cxform) is a library based on the algorithms in Hapgood, 1992. Python wrappers include:
[aics](https://aics.readthedocs.io/api.html), [python-magnetosphere](https://github.com/dpq/python-magnetosphere), and [bsd-conquerer](https://github.com/bsd-conqueror/cxform).

[Geopack-08 (Fortran)](https://ccmc.gsfc.nasa.gov/models/Tsyganenko%20Magnetic%20Field~T96/) is a library of utility functions related to magnetospheric magnetic field models. Python translations and wrappers include [Geopack (Python)](https://pypi.org/project/geopack/) is based on a hand translation of Tysganenko's Geopack (Fortran) to native Python and [PyGeopack](https://pypi.org/project/PyGeopack/), which wraps Geopack using Python `ctypes` and requires the user to provide a compiled Geopack shared object library or DLL.

`hxform` includes wrappers to [Tsyganenko's Geopack-08 library](https://ccmc.gsfc.nasa.gov/models/Tsyganenko%20Magnetic%20Field~T96/) and [cxform](https://github.com/edsantiago/cxform), both of which contain magnetospheric coordinate transformation functions. Copies of the source code are in [src](https://github.com/rweigel/hxform/tree/master/src)).

To wrap `Geopack-08`, Numpy's `f2py` is used; see [`Geopack-2008_dp_wrapper.for`](https://github.com/rweigel/hxform/tree/master/src/geopack-08).

To wrap `cxform` Python's `ctype` library is used; see [`cxform_wrapper.c`](https://github.com/rweigel/hxform/tree/master/src/cxform). Arrays are passed to the wrapper function, which loops over the array and calls the required functions on each iteration. This is much faster than looping over an array in Python and calling an external library function on each iteration.

To build the libraries manually, install `NumPy` and `meson-python` and

```bash
cd src/geopack-2008
# The following will generate many warnings
f2py -c Geopack-2008_dp_wrapper.for Geopack-2008_dp.for T96_01.for -m geopack_08_dp
cp *.so ../../hxform
```

```bash
cd src/cxform/
gcc -fPIC -shared -o cxform_wrapper.so cxform_wrapper.c cxform-manual.c cxform-auto.c
cp *.so ../../hxform
```

The following demos show how `hxform` calls the Geopack-08 and cxform wrappers. In general, you should not need to use these methods except for debugging.

```bash
python demo-native/geopack_08_dp_demo.py
python demo-native/cxform_demo.py
```