# Compare differences and timing

arr_size = 2
print_all = False

import os
import sys
import time as tm
import numpy as np
from datetime import datetime

import warnings
warnings.filterwarnings("ignore", message="Leapseconds.*")

from hxform import hxform as hx

from xprint import Xprint as Xp
xp = Xp() # Print to console and log file

X1 = (200/arr_size)*np.arange(0.005, arr_size, dtype=np.float64)
Y1 = (200/arr_size)*np.arange(0.005, arr_size, dtype=np.float64)
Z1 = (200/arr_size)*np.arange(0.005, arr_size, dtype=np.float64)

xforms = []
csys = ['GEI','GEO','GSE','GSM','MAG','SM']
#csys = ['GEO','GSE']

for c1 in csys:
    for c2 in csys:
        if c1 != c2:
            xforms.append('to'.join((c1, c2)))

XYZ = np.column_stack([X1,Y1,Z1])
time = (2000,1,10)
#time = np.full(XYZ.shape, time)

for xform in xforms:

    initial, final = xform.split('to')

    gp_start = tm.time()
    out_gp = hx.transform(XYZ, time, initial, final, lib='geopack_08_dp')
    gp_end = tm.time()

    sp_start = tm.time()
    out_sp = hx.transform(XYZ, time, initial, final, lib='spacepy')
    sp_end = tm.time()

    cx_start = tm.time()
    out_cx = hx.transform(XYZ, time, initial, final, lib='cxform')
    cx_end = tm.time()

    magdiff = np.linalg.norm(out_sp-out_gp)/np.linalg.norm(out_gp)
    maxdiff_sp_gp = np.max((out_sp - out_gp)/out_gp, axis=0)

    magdiff = np.linalg.norm(out_sp-out_gp)/np.linalg.norm(out_gp)
    maxdiff_cx_gp = np.max((out_cx - out_gp)/out_gp, axis=0)

    gp_time = gp_end - gp_start
    sp_time = sp_end - sp_start
    cx_time = cx_end - cx_start

    sp_gp_err = ""
    if np.any( np.abs(maxdiff_sp_gp) > 0.001):
        sp_gp_err = "!!!"
    cx_gp_err = ""
    if np.any( np.abs(maxdiff_cx_gp) > 0.001):
        cx_gp_err = "!!!"

    if print_all == False and (sp_gp_err == "" and cx_gp_err == ""):
        continue

    xp.xprint(60*"-")
    xp.xprint("Transform: {0:s}".format(xform))
    xp.xprint("Time: {0:d}-{1:02d}-{2:02d}".format(*time))

    xp.xprint("# of points: {0:d}".format(arr_size))
    xp.xprint("First input: [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[0], Y1[0], Z1[0]))
    xp.xprint("Last input:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[-1], Y1[-1], Z1[-1]))

    if True:

        xp.xprint("Results:")
        xp.xprint("    SpacePy")
        xp.xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_sp[0][0], out_sp[0][1], out_sp[0][2]))
        xp.xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_sp[-1][0], out_sp[-1][1], out_sp[-1][2]))
        xp.xprint("      Run time: {0:.2e} s".format(sp_time))
        xp.xprint("    Geopack-08 double precision")
        xp.xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_gp[0][0], out_gp[0][1], out_gp[0][2]))
        xp.xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_gp[-1][0], out_gp[-1][1], out_gp[-1 ][2]))
        xp.xprint("      Run time: {0:.2e} s".format(gp_time))
        xp.xprint("    cxform")
        xp.xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_cx[0][0], out_cx[0][1], out_cx[0][2]))
        xp.xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]"\
            .format(out_cx[-1][0], out_cx[-1][1], out_cx[-1 ][2]))
        xp.xprint("      Run time: {0:.2e} s".format(gp_time))

    #print(out_gp)
    #print(out_sp)
    #print(out_cx)

    xp.xprint("SpacePy/Geopack-08-dp time: {0:.1f}".format(sp_time/gp_time))
    xp.xprint("cxform/Geopack-08-dp time:  {0:.1f}".format(cx_time/gp_time))
    xp.xprint(sp_gp_err + "Max (SpacePy - Geopack-08-dp)/Geopack-08-dp: [{0:.2e}, {1:.2e}, {2:.2e}]"\
                .format(maxdiff_sp_gp[0], maxdiff_sp_gp[1], maxdiff_sp_gp[2]))
    xp.xprint(cx_gp_err + "Max (cxform  - Geopack-08-dp)/Geopack-08-dp: [{0:.2e}, {1:.2e}, {2:.2e}]"\
                .format(maxdiff_cx_gp[0], maxdiff_cx_gp[1], maxdiff_cx_gp[2]))
    xp.xprint(60*"-")
