# hxform

Heliophysics coordinate transforms in Python

# Overview

This package is a thin and fast wrapper to [Tsyganenko's Geopack-08 library](https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=Tsyganenko%20Magnetic%20Field), which contains magnetospheric coordinate tranformation functions. To wrap this library, Numpy's `f2py` is used; see `src/Geopack-2008_dp_wrapper.for`. Arrays are passed to the Fortran wrapper function, which loops over it and calls the required Geopack functions on each iteration. This is much faster than looping over an array in Python and calling a Geopack functions on each iteration.

`hxform` also contains a wrapper to [SpacePy's coordinate tranformation functions](https://spacepy.github.io/irbempy.html), which requires the installation of SpacePy. (SpacePy is not installed when `hxform` is installed due to issues encountered with SpacePy installation at the time of this release.)

Extensive testing and inter-comparison has been performed on coordinate transform calculations between this library, [SpacePy (Python)](https://spacepy.github.io/irbempy.html) (which wraps [IRBEM (Fortran)](https://sourceforge.net/projects/irbem/), which in turn wraps Tsyganenko's Geopack (Fortran)), and [SSCWeb's coordinate calculator](https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi) (which uses [CXFORM](https://spdf.gsfc.nasa.gov/pub/software/old/selected_software_from_nssdc/coordinate_transform/)).

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

Alternative approach:

```bash
# The following will generate many warnings
f2py -c Geopack-2008_dp_wrapper.for Geopack-2008_dp.for T96_01.for -m geopack_08_dp
cp *.so ../hxform
cd ..;
pip install -e .
```

# Tests, Comparisons, and Demos

See the files in the [test directory](https://github.com/rweigel/hxform/tree/master/test). The result form executing the files is stored in a `.log` file.

```
python test/hxform_test1.py
python test/hxform_test2.py
python test/geopack_08_dp_wrapper_demo.py
```

Demo of bug in SpacePy

```
python test/spacepy_demo.py
```

# Related Code

* [SpacePy (Python)](https://spacepy.github.io/irbempy.html) wraps [IRBEM (Fortran)](https://sourceforge.net/projects/irbem/) which wraps Tsyganenko's GEOPACK (Fortran) using Python ctypes (check this). 
* [CXFORM (C)](https://spdf.gsfc.nasa.gov/pub/software/old/selected_software_from_nssdc/coordinate_transform/) is a library used by [SSCWeb](https://sscweb.gsfc.nasa.gov/) for coordinate tranformations. The library is based on the algorithms in Hapgood, 1992. Python wrappers include:
  * https://aics.readthedocs.io/api.html
  * https://github.com/dpq/python-magnetosphere
* [Geopack (Python)](https://pypi.org/project/geopack/) is based on a hand translation of Tysganenko's Geopack (Fortran) to native Python.
* [PyGeopack](https://pypi.org/project/PyGeopack/) wraps Geopack using Python `ctypes` and requires the user to provide a compiled Geopack shared object library or DLL.

# References

* [Laundal and Richmond, 2016, Magnetic Coordinate Systems](https://arxiv.org/ct?url=https%3A%2F%2Fdx.doi.org%2F10.1007%2Fs11214-016-0275-y&v=34afcdf3)
* [Hapgood, 1992, Space physics coordinate transformations: A user guide](https://doi.org/10.1016/0032-0633(92)90012-D)
* [Kivelson and Russell, 1995; Appendix 3](https://books.google.com/books/about/Introduction_to_Space_Physics.html?id=qWHSqXGfsfQC)
* [Russell, 1971, Geophysical Coordinate Tranformations](http://jsoc.stanford.edu/~jsoc/keywords/Chris_Russel/Geophysical%20Coordinate%20Transformations.htm)
* [SPENVIS help on coordinate transformations](https://www.spenvis.oma.be/help/background/coortran/coortran.html); includes [an animation](https://www.spenvis.oma.be/help/background/coortran/anim.html) and additional references at the bottom of the page.

