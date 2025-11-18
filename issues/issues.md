# Reports

* https://github.com/spacepy/spacepy/issues/534
* https://github.com/sunpy/sunpy/issues/8189
* https://github.com/sunpy/sunpy/issues/8188
* https://github.com/sunpy/sunpy/issues/8406
* https://github.com/spedas/pyspedas/issues/1207
* https://github.com/spedas/pyspedas/issues/1273
* https://github.com/spedas/pyspedas/issues/1274

# Geopack

## f2py compilation error with EXNAME and INNAME
Compiling with f2py caused errors with properly identifying EXNAME and INNAME which are function variable names passed in as arguments into subroutines TRACE\_08, RHAND\_08, and STEP\_08. The error did not allow compilation to complete.

**solution** 
- The variables EXNAME and INNAME were modified to be strings. Then depending on the string value the correct function name would be called.


