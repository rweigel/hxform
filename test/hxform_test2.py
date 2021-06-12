arr_size = 1000

import os
import sys
import time as tm
import numpy as np
from datetime import datetime

from hxform import hxform as hx

# Print to console and logfile
import os
logfile = os.path.realpath(__file__)[0:-2] + "log"
if os.path.exists(logfile):
    with open(logfile, "w") as f: pass
def xprint(msg):
    import os; print(msg);
    with open(logfile, "a") as f: f.write(str(msg) + "\n")

X1 = (200/arr_size)*np.arange(0, arr_size, dtype=np.float64)
Y1 = (200/arr_size)*np.arange(0, arr_size, dtype=np.float64)
Z1 = (200/arr_size)*np.arange(0, arr_size, dtype=np.float64)

xforms = []
#csys = ['GSM','GSE','GSE','MAG','GEO','GEI']
csys = ['GSM','GSE']
for c1 in csys:
    for c2 in csys:
        if c1 != c2:
            xforms.append('to'.join((c1, c2)))

XYZ = np.column_stack([X1,Y1,Z1])
time = (1997,1,1)
#time = np.full(XYZ.shape, time)

for xform in xforms:

    initial, final = xform.split('to')

    gp_start = tm.time()
    out_gp = hx.transform(XYZ, time, initial, final, lib='geopack_08_dp')
    gp_end = tm.time()

    sp_start = tm.time()
    out_sp = hx.transform(XYZ, time, initial, final, lib='spacepy')
    sp_end = tm.time()

    magdiff = np.linalg.norm(out_gp-out_sp)/np.linalg.norm(out_sp)
    maxdiff = np.max(out_gp - out_sp, axis=0)

    gp_time = gp_end - gp_start
    sp_time = sp_end - sp_start

    xprint(60*"-")
    xprint("transform:   {0:s}".format(xform))
    xprint("# of points: {0:d}".format(arr_size))
    xprint("First input: [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[0], Y1[0], Z1[0]))
    xprint("Last input:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[-1], Y1[-1], Z1[-1]))
    #xprint("Time:        {0:s}".format(time_str))

    xprint("Results:")
    xprint("    SpacePy")
    xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_sp[0][0], out_sp[0][1], out_sp[0][2]))
    xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_sp[-1][0], out_sp[-1][1], out_sp[-1][2]))
    xprint("      Run time: {0:.2e} s".format(sp_time))
    xprint("    Geopack-08 double precision")
    xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_gp[0][0], out_gp[0][1], out_gp[0][2]))
    xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_gp[-1][0], out_gp[-1][1], out_gp[-1 ][2]))
    xprint("      Run time: {0:.2e} s".format(gp_time))

    xprint("SpacePy/Geopack-08-dp time: {0:.1f}".format(sp_time/gp_time))
    xprint("Max output difference: [{0:.2e}, {1:.2e}, {2:.2e}]".format(maxdiff[0], maxdiff[1], maxdiff[2]))
    xprint(60*"-")
