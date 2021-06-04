# Errors and solutions encountered

## spacepy


### GSM to SM produces nan values

transforming GSM to SM produces an error. Other transformations including SM to GSM do not produce errors. This error did not occur in the past.  

```python 

import spacepy.coordinates as sc
from spacepy.time import Ticktock
import spacepy
import sys

print("spacepy version: ",spacepy.__version__)
print("python version: ", sys.version)

cvals = sc.Coords([[1,2,4],[1,2,2]], 'GSM', 'car')
cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticks
newcoord = cvals.convert('SM', 'car')
print(newcoord)

```

The following was the output. 

```
spacepy version:  0.2.2
python version:  3.8.10 | packaged by conda-forge | (default, May 11 2021, 06:39:48) 
[Clang 11.1.0 ]
Coords( [[nan, nan, nan], [nan, nan, nan]] , 'SM', 'car')
```

**no solution yet** 
- tried reinstalling spacepy inside of a fresh anaconda environment. 



## Geopack

### f2py compilation error with EXNAME and INNAME
Compiling with f2py caused errors with properly identifying EXNAME and INNAME which are function variable names passed in as arguments into subroutines TRACE\_08, RHAND\_08, and STEP\_08. The error did not allow compilation to complete. 

**solution** 
- The variables EXNAME and INNAME were modified to be strings. Then depending on the string value the correct function name would be called. 


