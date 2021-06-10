# Demo of calling Fortran library geopack_08_dp directly instead of through hxform.

import numpy as np
import hxform.geopack_08_dp as geopack_08_dp
import hxform.hxform as hx

def xprint(msg):
    # Print to console and logfile
    import os
    print(msg);
    logfile = os.path.realpath(__file__)[0:-2] + "log"
    if not os.path.exists(logfile):
        with open(logfile, "w") as f: pass
    with open(logfile, "a") as f: f.write(str(msg) + "\n")

time = (1997,1,1)
input = [1., 0., 0.]
# Time format for hxform is
# year, month, day, [hours, [minutes, [seconds]]]
# only year, month, day are required.
output = hx.GSMtoSM(input, time, lib='geopack_08_dp')
xprint(input)
xprint(output)

# Time format for geopack_08_dp is
# year, doy, hour, minute, second
# all values are required.
dtime = np.array((1997,1,1,0,0), dtype=np.int32)
input = [1., 0., 0.]
output = np.column_stack(geopack_08_dp.transform(input[0], input[1], input[2], "GSMtoSM", dtime))

xprint(input)
xprint(output)

input = [[1., 0., 0.],[1., 0., 0.]]
output = hx.GSMtoSM(input, time, lib='geopack_08_dp')
xprint(input)
xprint(output)

dtime = np.array((1997,1,1,0,0), dtype=np.int32)
input = [[1., 1., 1.],[0., 0., 0.],[0., 0., 0.]]
output = np.column_stack(geopack_08_dp.transform(input[0], input[1], input[2], "GSMtoSM", dtime))

xprint(input)
xprint(output)

