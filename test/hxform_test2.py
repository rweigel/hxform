arr_size = 10000

import os
import sys
import time as tm
import numpy as np
from datetime import datetime

import hxform.hxform as hx

def xprint(msg):
    # Print to console and logfile
    import os
    print(msg);
    logfile = os.path.realpath(__file__)[0:-2] + "log"
    if not os.path.exists(logfile):
        with open(logfile, "w") as f: pass
    with open(logfile, "a") as f: f.write(msg + "\n")

# Tsyganenko uses Geocentric-Solar Wind GSW instead of GSM. GSW has the positive
# x-axis pointing antiparallel to the solar wind. He chose this coordinate system
# in order to simplify his code. GSW becomes identical to GSM when the solar Wind
# velocity travels only in the negative x-axis direction, Tsyganenko chose
#  <-400,0,0> SMGSW_08

X1 = (200/arr_size)*np.arange(0, arr_size, dtype=np.float64)
Y1 = (200/arr_size)*np.arange(0, arr_size, dtype=np.float64)
Z1 = (200/arr_size)*np.arange(0, arr_size, dtype=np.float64)

trans = []
csys = ['GSM','SM']#,'GSE','MAG','GEO','GEI']
for c1 in csys:
    for c2 in csys:
        if c1 != c2:
            trans.append('to'.join((c1,c2)))

XYZ = np.column_stack([X1,Y1,Z1])
time = (1997,1,1)
time_str = datetime(time[0], time[1], time[2]).isoformat()

for t in trans:

    initial, final = t.split('to')

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
    xprint("transform:   {0:s}".format(t))
    xprint("# of points: {0:d}".format(arr_size))
    xprint("First input: [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[0], Y1[0], Z1[0]))
    xprint("Last input:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[-1], Y1[-1], Z1[-1]))
    xprint("Time:        {0:s}".format(time_str))

    xprint("Results:")
    xprint("    SpacePy")
    xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_sp[0][0], out_sp[0][1], out_sp[0][2]))
    xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_sp[-1][0], out_sp[-1][1], out_sp[-1][2]))
    xprint("      Run time: {0:.2e} s".format(sp_time))
    xprint("    hxform")
    xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_gp[0][0], out_gp[0][1], out_gp[0][2]))
    xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(out_gp[-1][0], out_gp[-1][1], out_gp[-1 ][2]))
    xprint("      Run time: {0:.2e} s".format(gp_time))

    xprint("SpacePy/hxform time:   {}".format(sp_time/gp_time))
    xprint("Max output difference: [{0:.2e}, {1:.2e}, {2:.2e}]".format(maxdiff[0], maxdiff[1], maxdiff[2]))
    xprint(60*"-")
