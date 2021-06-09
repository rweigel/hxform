import numpy as np
from datetime import datetime

import hxform.cxtransform as cx
import hxform.hxform as hx

def xprint(msg):
    # Print to console and logfile
    import os
    print(msg);
    logfile = os.path.realpath(__file__)[0:-2] + "log"
    if not os.path.exists(logfile):
        with open(logfile, "w") as f: pass
    with open(logfile, "a") as f: f.write(msg + "\n")

# [year, month, day, hours, minutes, seconds, x, y, z, 'car']
# or
# [year, month, day, hours, minutes, seconds, r, long in degrees, latitude in degrees, 'sph']

data = [ 
        [2003.0, 11.0, 20.0, 7.0,   0.0, 0.0, 1.,   0. ,  0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,   0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,  68.5, 50.0, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 43.0, 0.0, 1., 166.0, 52.5, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 47.0, 0.0, 1., 166.0, 52.5, 'sph']
    ]

# spherical input/output option note implemented in hxtranform with lib=geopack_08_dp
# So for now, only use these two checks.
data = [ 
        [2003.0, 11.0, 20.0, 7.0,   0.0, 0.0, 1.,   0. ,  0.,  'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,   0.,  'car']
    ]

# Results from
# https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi
# Can only manually enter spherical coordinates for MAG

# GSM [x, y, z, MLT]
sscweb = [
            [-0.68, -0.64, -0.36,  2.64457],
            [ 0.94,  0.30,  0.16, 13.16680],
            [-0.09,  0.64,  0.76, 17.73347],
            [-0.72, -0.18,  0.67,  1.13969],
            [-0.72, -0.19,  0.67,  1.20329]
        ]


#lib = 'spacepy'
lib = 'geopack_08_dp'

k = 0
for d in data:
    print('------------------------------------------------------------')
    syst = d[9] # sph or cart
    time = d[0:6]
    if syst == 'car':
        pos = d[6:9]
        v_sp = hx.MAGtoGSM(pos, time, ctype_in='car', ctype_out='car', lib='spacepy')
        MLT_sp = hx.MAGtoMLT(pos, [time], csys='car')
        v_gp = hx.MAGtoGSM(pos, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp')
        MLT_gp = hx.MAGtoMLT(pos, [time], csys='car')
    elif syst == 'sph':
        r, MLON, MLAT = d[6:9]
        v = hx.MAGtoGSM([r, MLAT, MLON], time, ctype_in='sph', ctype_out='car', lib=lib)
        v = cx.MAGtoGSM([r, MLAT, MLON], [time], 'sph', 'car')
        MLT = cx.MAGtoMLT(MLON, [time])
    else:
        print('INVALID COORDINATE TYPE. Use "car" or "sph"')
    UT = time[3] + time[4]/60.+ time[5]/3600.

    # NOTE: if MLT has more then one element, this will be longer than 4 and only 
    # the first element of MLT will show up in the print statements that follow
    spacepy = np.append(v_sp, MLT_sp)
    geopack = np.append(v_gp, MLT_gp)
    #print(spacepy.shape)

    time_str = datetime(int(time[0]), int(time[1]), int(time[2]), int(time[3]), int(time[4]), int(time[5])).isoformat()

    xprint('Input:')
    xprint('   ' + time_str)
    if syst == 'car':
        xprint('   MAG           x     y     z')
    else:
        xprint('   MAG           r     mlon    mlat')
    xprint('                {0:.2f}  {1:.2f}  {2:.2f}'.format(d[6], d[7], d[8]))
    xprint('Output:')
    xprint('   GSM              x           y           z          MLT')
    xprint('   SpacePy:    {0:.8f}  {1:.8f}  {2:.8f}  {3:11.8f}' \
        .format(spacepy[0], spacepy[1], spacepy[2], spacepy[3]))
    xprint('   Geopack:    {0:.8f}  {1:.8f}  {2:.8f}  {3:11.8f}' \
        .format(geopack[0], geopack[1], geopack[2], geopack[3]))
    xprint('   SSCWeb:     {0:.8f}  {1:.8f}  {2:.8f}  {3:11.8f}' \
        .format(sscweb[k][0], sscweb[k][1], sscweb[k][2], sscweb[k][3]))
    #print('   Diff:       {0:.3f} {1:.3f} {2:.3f} {3:.3f}' \
    #    .format(sscweb[k][0] - spacepy[0],sscweb[k][1] - spacepy[1],
    #            sscweb[k][2] - spacepy[2],sscweb[k][3] - spacepy[3]))
    k = k + 1
