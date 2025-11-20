# hxform

Heliophysics coordinate transforms in Python using various libraries, Python packages, and web services using a common interface.

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

See the [`demo/`](https://github.com/rweigel/hxform/tree/main/demo) directory for examples.

# Overview

The list of supported libraries and frames is in [`info.py`](https://github.com/rweigel/hxform/blob/main/hxform/info.py).

This package includes wrappers to [Tsyganenko's Geopack-08 library](https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=Tsyganenko%20Magnetic%20Field) and [cxform](https://github.com/edsantiago/cxform), both of which contain magnetospheric coordinate transformation functions.

To wrap `Geopack-08`, Numpy's `f2py` is used; see `src/Geopack-2008_dp_wrapper.for`. To wrap `cxform` Python's `ctype` library is used; see `src/cxform_wrapper.c`. Arrays are passed to the wrapper function, which loops over the array and calls the required functions on each iteration. Looping in Fortran or c is much faster than looping over an array in Python and calling an external library function on each iteration.

# Install

Requires Anaconda and a Fortran compiler.

Tested with Python 3.9.12, 3.10.14, and 3.11.9. (Install fails in Python 3.12 due to `numpy.distutils` and `distutils` deprecation.)

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

See the files in the [test](https://github.com/rweigel/hxform/tree/master/test) directory. The result from executing the files is stored in a `.log` file.

To run, execute

```bash
pip install pytest
pytest # Executes files in ./test
```

# Development

To build the libraries for testing, use

```bash
cd src/geopack-2008
# The following will generate many warnings
f2py -c Geopack-2008_dp_wrapper.for \
 Geopack-2008_dp.for \
        T96_01.for \
 -m geopack_08_dp
cp *.so ../../hxform
```

```bash
cd src/cxform/
gcc -fPIC -shared -o cxform_wrapper.so cxform_wrapper.c cxform-manual.c cxform-auto.c
cp *.so ../../hxform
```

The following demos show how `hxform` calls Geopack-08 and cxform. In general, you should not need to use these methods except for debugging.

```bash
python demo-native/geopack_08_dp_wrapper_demo.py
python demo-native/cxform_demo.py
```

Examples of using the native interface for other Python packages are in the [`demo-native`](https://github.com/rweigel/hxform/tree/main/demo-native) directory.

# Related Code

* [cxform (C)](https://github.com/edsantiago/cxform) is a library based on the algorithms in Hapgood, 1992. Python wrappers include:
  * https://aics.readthedocs.io/api.html
  * https://github.com/dpq/python-magnetosphere
  * https://github.com/bsd-conqueror/cxform
* [Geopack-08 (Fortran)](https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=Tsyganenko%20Magnetic%20Field) is a library of utility functions related to magnetospheric magnetic field models. Python translations and wrappers include
  * [Geopack (Python)](https://pypi.org/project/geopack/) is based on a hand translation of Tysganenko's Geopack (Fortran) to native Python. 
  * [PyGeopack](https://pypi.org/project/PyGeopack/) wraps Geopack using Python `ctypes` and requires the user to provide a compiled Geopack shared object library or DLL.
  * [SpacePy (Python)](https://spacepy.github.io/irbempy.html), which wraps [IRBEM (Fortran)](https://sourceforge.net/projects/irbem/).

# References

* [Laundal and Richmond, 2016, Magnetic Coordinate Systems](https://arxiv.org/ct?url=https%3A%2F%2Fdx.doi.org%2F10.1007%2Fs11214-016-0275-y&v=34afcdf3)
* [Hapgood, 1992, Space physics coordinate transformations: A user guide](https://doi.org/10.1016/0032-0633(92)90012-D)
* [Kivelson and Russell, 1995; Appendix 3](https://books.google.com/books/about/Introduction_to_Space_Physics.html?id=qWHSqXGfsfQC)
* [Russell, 1971, Geophysical Coordinate Transformations](http://jsoc.stanford.edu/~jsoc/keywords/Chris_Russel/Geophysical%20Coordinate%20Transformations.htm)
* [SPENVIS help on coordinate transformations](https://www.spenvis.oma.be/help/background/coortran/coortran.html); includes [an animation](https://www.spenvis.oma.be/help/background/coortran/anim.html) and additional references at the bottom of the page.
