# Demo of calling Fortran library geopack_08_dp directly instead of through hxform.

import numpy as np
import hxform.geopack_08_dp as geopack_08_dp
import hxform.hxform as hx

from hxform import xprint # print to console and geopack_08_dp_demo.log

# Time format for geopack_08_dp is
# year, doy, hour, minute, second
# all values are required.
dtime = np.array([[1997, 1, 1, 0, 0]], dtype=np.int32)
input = np.array([[1., 0., 0.]], dtype=np.float64)
output = geopack_08_dp.transform(input, "GSMtoGSE", dtime, 1)
