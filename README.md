# hxform

Heliophysics coordinate transforms in Python

# Overview

This package is a thin and fast wrapper to [Tsyganenko's Geopack-08 library](https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=Tsyganenko%20Magnetic%20Field) and [cxform](https://github.com/edsantiago/cxform), both of which contains magnetospheric coordinate tranformation functions. 

To wrap ``Geopack-08`, Numpy's `f2py` is used; see `src/Geopack-2008_dp_wrapper.for`. To wrap `csform` Python's `ctype` library is used; see `src/cxform_wrapper.for`. For both wrappers, arrays are passed to the  wrapper function, which loops over it and calls the required functions on each iteration. This is much faster than looping over an array in Python and calling an external library function on each iteration.

`hxform` also contains a wrapper to [SpacePy's coordinate tranformation functions](https://spacepy.github.io/irbempy.html), which requires the installation of SpacePy. (SpacePy is not installed when `hxform` is installed due to issues encountered with SpacePy installation at the time of this release.)

# Install

```bash  
git clone https://github.com/rweigel/hxform
cd hxform;
pip install -e .
```

To test if installation was successful, execute

```
python hxform_demo.py
```

See also the files in the [demo](https://github.com/rweigel/hxform/tree/master/demo) directory. The result form executing the files is stored in a `.log` file.

# Tests

See the files in the [test](https://github.com/rweigel/hxform/tree/master/test). The result form executing the files is stored in a `.log` file.

To run, execute

```bash
python test/hxform_test1.py
python test/hxform_test2.py
```

# Development

To build the libraries for testing, use

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

The following demos show how [`hxform.py`](https://github.com/rweigel/hxform/blob/main/hxform/hxform.py) calls Geopack-08, cxform, and SpacePy. In general, you should not need to use these methods except for debugging.

```bash
python demo/geopack_08_dp_wrapper_demo.py
python demo/cxform_demo.py
python demo/spacepy_demo.py
```
# Related Code

* [cxform (C)](https://github.com/edsantiago/cxform) is a library is based on the algorithms in Hapgood, 1992. Python wrappers include:
  * https://aics.readthedocs.io/api.html
  * https://github.com/dpq/python-magnetosphere
  * https://github.com/bsd-conqueror/cxform
* [Geopack-08 (Fortran)](https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=Tsyganenko%20Magnetic%20Field) is a library of utility functions related to magnetospheric magentic field models. Python translations and wrappers include
  * [Geopack (Python)](https://pypi.org/project/geopack/) is based on a hand translation of Tysganenko's Geopack (Fortran) to native Python. 
  * [PyGeopack](https://pypi.org/project/PyGeopack/) wraps Geopack using Python `ctypes` and requires the user to provide a compiled Geopack shared object library or DLL.
  * [SpacePy (Python)](https://spacepy.github.io/irbempy.html) wraps [IRBEM (Fortran)](https://sourceforge.net/projects/irbem/) which wraps Tsyganenko's GEOPACK (Fortran) using Python ctypes (check this). 

# References

* [Laundal and Richmond, 2016, Magnetic Coordinate Systems](https://arxiv.org/ct?url=https%3A%2F%2Fdx.doi.org%2F10.1007%2Fs11214-016-0275-y&v=34afcdf3)
* [Hapgood, 1992, Space physics coordinate transformations: A user guide](https://doi.org/10.1016/0032-0633(92)90012-D)
* [Kivelson and Russell, 1995; Appendix 3](https://books.google.com/books/about/Introduction_to_Space_Physics.html?id=qWHSqXGfsfQC)
* [Russell, 1971, Geophysical Coordinate Tranformations](http://jsoc.stanford.edu/~jsoc/keywords/Chris_Russel/Geophysical%20Coordinate%20Transformations.htm)
* [SPENVIS help on coordinate transformations](https://www.spenvis.oma.be/help/background/coortran/coortran.html); includes [an animation](https://www.spenvis.oma.be/help/background/coortran/anim.html) and additional references at the bottom of the page.

