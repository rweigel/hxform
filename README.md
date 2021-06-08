# hxtransform

Heliophysics coordinate transforms in Python

This package is a thin wrapper to Tsyganenko's GEOPACK library, which contains magnetospheric coordinate tranformation functions. Although SpacePy also wraps this library (via another wrapper library), its wrapper is not implemented for optimial efficiency - when an array of coordinates are to be transformed, the coordinates are looped over in Python and the Fortran function is called on each iteration. In contrast, this library performs the looping in Fortran, which leads to 100x or larger speed-ups, which was needed for the use-case of transforming the output of magnetospheric MHD simulations.

# Related Code

* [SpacePy (Python)](https://spacepy.github.io/irbempy.html) wraps [IRBEM (Fortran)](https://sourceforge.net/projects/irbem/) which wraps Tsyganenko's GEOPACK (Fortran).
* [CXFORM (C)](https://spdf.gsfc.nasa.gov/pub/software/old/selected_software_from_nssdc/coordinate_transform/) is a library use by [SSCWeb](https://sscweb.gsfc.nasa.gov/) for coordinate tranformations. The library is based on the algorithms in Hapgood, 1992.

# References

* [Hapgood, 1992, Space physics coordinate transformations: A user guide](https://doi.org/10.1016/0032-0633(92)90012-D)
* [Kivelson and Russell, 1995; Appendix 3](https://books.google.com/books/about/Introduction_to_Space_Physics.html?id=qWHSqXGfsfQC)
* [Russell, 1971, Geophysical Coordinate Tranformations](http://jsoc.stanford.edu/~jsoc/keywords/Chris_Russel/Geophysical%20Coordinate%20Transformations.htm)
* [SPENVIS help on coordinate transformations](https://www.spenvis.oma.be/help/background/coortran/coortran.html]; includes (an animation)[https://www.spenvis.oma.be/help/background/coortran/coortran.html) and additional references at the bottom of the page.

# Installation

To compile hxtransform python package, execute
```bash  
f2py -c Geopack_08_dp.for T96_01.for -m hxtransform
```

To test if compilation was successful, execute
``` python -c "import hxtransform as hx; help(hx)"```

Which should produce something similar to 
```
<module 'hxtransform' from '/Users/user/hxtransform/hxtransform.cpython-38-darwin.so'> 
``` 
