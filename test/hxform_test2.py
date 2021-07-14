arr_size = 10000

import os
import sys
import time as tm
import numpy as np
from datetime import datetime

import warnings
warnings.filterwarnings("ignore", message="Leapseconds.*")

from hxform import hxform as hx

from hxform.xprint import Xprint as Xp
xp = Xp() # Print to console and log file

X1 = np.linspace(1., 200., arr_size, dtype=np.float64)
Y1 = np.linspace(1., 200., arr_size, dtype=np.float64)
Z1 = np.linspace(1., 200., arr_size, dtype=np.float64)

xforms = []
csys = ['GEI','GEO','GSE','GSM','MAG','SM']
#csys = ['GSM','GSE']

for c1 in csys:
    for c2 in csys:
        if c1 != c2:
            xforms.append('to'.join((c1, c2)))

XYZ = np.column_stack([X1,Y1,Z1])
time = (1997,1,10)
#time = np.full(XYZ.shape, time)

for xform in xforms:

    initial, final = xform.split('to')

    gp_start = tm.time()
    out_gp = hx.transform(XYZ, time, initial, final, lib='geopack_08_dp')
    gp_end = tm.time()

    sp_start = tm.time()
    out_sp = hx.transform(XYZ, time, initial, final, lib='spacepy')
    sp_end = tm.time()

    sp2_start = tm.time()
    out_sp2 = hx.transform(XYZ, time, initial, final, lib='spacepy-irbem')
    sp2_end = tm.time()

    cx_start = tm.time()
    out_cx = hx.transform(XYZ, time, initial, final, lib='cxform')
    cx_end = tm.time()

    maxdiff_sp_gp = np.max((out_sp - out_gp)/np.linalg.norm(out_gp), axis=0)
    maxdiff_sp2_gp = np.max((out_sp2 - out_gp)/np.linalg.norm(out_gp), axis=0)
    maxdiff_cx_gp = np.max((out_cx - out_gp)/np.linalg.norm(out_gp), axis=0)

    gp_time = gp_end - gp_start
    sp_time = sp_end - sp_start
    sp2_time = sp2_end - sp2_start
    cx_time = cx_end - cx_start

    xp.xprint(60*"-")
    xp.xprint("Transform: {0:s}".format(xform))
    xp.xprint("Time: {0:d}-{1:02d}-{2:02d}".format(*time))

    xp.xprint("# of points: {0:d}".format(arr_size))
    xp.xprint("First input: [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[0], Y1[0], Z1[0]))
    xp.xprint("Last input:  [{0:.8f}, {1:.8f}, {2:.8f}]".format(X1[-1], Y1[-1], Z1[-1]))

    if False:

        xp.xprint("Results:")
        xp.xprint("    SpacePy-New")
        xp.xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_sp[0][0], out_sp[0][1], out_sp[0][2]))
        xp.xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_sp[-1][0], out_sp[-1][1], out_sp[-1][2]))
        xp.xprint("    SpacePy-IRBEM")
        xp.xprint("      First output: [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_sp2[0][0], out_sp2[0][1], out_sp2[0][2]))
        xp.xprint("      Last output:  [{0:.8f}, {1:.8f}, {2:.8f}]"\
                .format(out_sp2[-1][0], out_sp2[-1][1], out_sp2[-1][2]))
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

    large_diff = False
    sp_gp_err = ""
    if np.any( np.abs(maxdiff_sp_gp) > 1e-5):
        sp_gp_err = "*"
        large_diff = True

    sp2_gp_err = ""
    if np.any( np.abs(maxdiff_sp2_gp) > 1e-5):
        sp2_gp_err = "*"
        large_diff = True

    cx_gp_err = ""
    if np.any( np.abs(maxdiff_cx_gp) > 1e-5):
        cx_gp_err = "*"
        large_diff = True

    if True:
        xp.xprint("A - B Relative Diff = max( (A - B)/|B| ); * = a relative diff > 1e-5")
        xp.xprint("SpacePy-New/Geopack-08-dp time:   {0:.1f}".format(sp_time/gp_time))
        xp.xprint("SpacePy-IRBEM/Geopack-08-dp time: {0:.1f}".format(sp2_time/gp_time))
        xp.xprint("cxform/Geopack-08-dp time:        {0:.1f}".format(cx_time/gp_time))
        xp.xprint("Max SpacePy New - Geopack-08-dp Relative Diff:   [{0:+.2e}, {1:+.2e}, {2:+.2e}]{3:s}"\
                    .format(maxdiff_sp_gp[0], maxdiff_sp_gp[1], maxdiff_sp_gp[2], sp_gp_err))
        xp.xprint("Max SpacePy IRBEM - Geopack-08-dp Relative Diff: [{0:+.2e}, {1:+.2e}, {2:+.2e}]{3:s}"\
                    .format(maxdiff_sp2_gp[0], maxdiff_sp2_gp[1], maxdiff_sp2_gp[2], sp2_gp_err))
        xp.xprint("Max cxform - Geopack-08-dp Relative Diff:        [{0:+.2e}, {1:+.2e}, {2:+.2e}]{3:s}"\
                    .format(maxdiff_cx_gp[0], maxdiff_cx_gp[1], maxdiff_cx_gp[2], cx_gp_err))
        xp.xprint(60*"-")
