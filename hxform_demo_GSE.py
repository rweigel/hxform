# time and input can be NumPy arrays
import hxform
from hxform import xprint

time = hxform.timelib.iso2ints('2012-06-30T00:00:00Z')
input = [1., 1., 1.]
initial = 'GEO'
final1 = 'GSE'
final2 = 'GSE_Z_PRIMARY'

lib = 'spiceypy2'

output1 = hxform.transform(input, time, initial, final1, lib=lib)
output2 = hxform.transform(input, time, initial, final2, lib=lib)

xprint("")
xprint(f"Time:   {time}")
xprint(f"Input:  {input}")
xprint(f"Transform: {initial} => {final1}")
xprint(f"Output: {output1}")
xprint(f"Transform: {initial} => {final2}")
xprint(f"Output: {output2}")

import numpy as np
mag1 = np.linalg.norm(output1)
mag2 = np.linalg.norm(output2)
angle = (180/np.pi)*np.arccos(np.dot(output1, output2)/(mag1*mag2))
print(f"angle = {angle:.2e} degrees")
print(f"angle = {angle/3600:.2e} arcseconds")
