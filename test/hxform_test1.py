import numpy as np
from datetime import datetime

from hxform import hxform as hx

from hxform.xprint import Xprint as Xp
xp = Xp() # Print to console and log file

# [year, month, day, hours, minutes, seconds, x, y, z, 'car']
# or
# [year, month, day, hours, minutes, seconds, r, latitude in degrees, long in degrees, 'sph']

data = [ 
        [2003.0, 11.0, 20.0, 7.0,   0.0, 0.0, 1.,   0.,    0., 'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1.,   0.,    0., 'car'],
        [2003.0, 11.0, 20.0, 17.0, 46.0, 0.0, 1., 50.0,  68.5, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 43.0, 0.0, 1., 52.5, 166.0, 'sph'],
        [2003.0, 11.0, 20.0, 18.0, 47.0, 0.0, 1., 52.5, 166.0, 'sph']
    ]

# Results from
# https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi
# Can only manually enter spherical coordinates for MAG

# GSM [x, y, z, MLT]
# For the first point, enter
# Time = 2003 324 17:46:00
# X = 1.    Lat: 0.
# Y = 0.    Lon: 0.
# Z = 0.    R: 1.
# Then click "GM" button (GM = MAG coordinate system)
# Output is
"""
TIME: 2003 324  7.00000 
              Input:  GM                                                        

 Radial distance     1.00 
          Lat    Long      X       Y       Z   hh.hhhhh 
 GEI    -10.29   92.13   -0.04    0.98   -0.18 
 J2000  -10.29   92.08   -0.04    0.98   -0.18 
 GEO    -10.29  -71.76    0.31   -0.93   -0.18  2.21613 
 GM       0.00    0.00    1.00    0.00    0.00  2.64457 
 GSE    -33.71 -145.02   -0.68   -0.48   -0.56  2.33229 
 GSM    -20.97 -136.87   -0.68   -0.64   -0.36 
 SM      -0.00 -140.33   -0.77   -0.64   -0.00  2.64457 
"""
# The hh.hhhh on the GM line is MLT
# The X, Y, Z values for GSM are MAG X,Y,Z = 1,0,0 converted to GSM

sscweb = [
            [-0.68, -0.64, -0.36,  2.64457],
            [ 0.94,  0.30,  0.16, 13.16680],
            [-0.09,  0.64,  0.76, 17.73347],
            [-0.72, -0.18,  0.67,  1.13969],
            [-0.72, -0.19,  0.67,  1.20329]
        ]

k = 0
for d in data:
    print('------------------------------------------------------------')
    syst = d[9] # sph or cart
    time = d[0:6]
    if syst == 'car':
        pos = d[6:9]
        v_sp = hx.MAGtoGSM(pos, time, ctype_in='car', ctype_out='car', lib='spacepy')
        MLT_sp = hx.MAGtoMLT(pos, time, csys='car', lib='spacepy')
        v_gp = hx.MAGtoGSM(pos, time, ctype_in='car', ctype_out='car', lib='geopack_08_dp')
        MLT_gp = hx.MAGtoMLT(pos, time, csys='car', lib='geopack_08_dp')
        v_cx = hx.MAGtoGSM(pos, time, ctype_in='car', ctype_out='car', lib='cxform')
        MLT_cx = hx.MAGtoMLT(pos, time, csys='car', lib='cxform')
    elif syst == 'sph':
        r, mlat, mlong = d[6:9]
        v_sp = hx.MAGtoGSM([r, mlat, mlong], time, ctype_in='sph', ctype_out='car', lib='spacepy')
        MLT_sp = hx.MAGtoMLT(mlong, time, csys='sph', lib='spacepy')
        v_gp = hx.MAGtoGSM([r, mlat, mlong], time, ctype_in='sph', ctype_out='car', lib='geopack_08_dp')
        MLT_gp = hx.MAGtoMLT(mlong, time, csys='sph', lib='geopack_08_dp')
        v_cx = hx.MAGtoGSM([r, mlat, mlong], time, ctype_in='sph', ctype_out='car', lib='cxform')
        MLT_cx = hx.MAGtoMLT(mlong, time, csys='sph', lib='cxform')
    else:
        print('INVALID COORDINATE TYPE. Use "car" or "sph"')
    UT = time[3] + time[4]/60.+ time[5]/3600.

    # NOTE: if MLT has more then one element, this will be longer than 4 and only 
    # the first element of MLT will show up in the print statements that follow
    spacepy = np.append(v_sp, MLT_sp)
    geopack = np.append(v_gp, MLT_gp)
    cxform = np.append(v_cx, MLT_cx)
    #print(spacepy.shape)

    time_str = datetime(int(time[0]), int(time[1]), int(time[2]), int(time[3]), int(time[4]), int(time[5])).isoformat()

    xp.xprint('Input:')
    xp.xprint('   ' + time_str)
    if syst == 'car':
        xp.xprint('   MAG           x     y     z')
    else:
        xp.xprint('   MAG           r     mlat   mlon')
    xp.xprint('                {0:.2f}  {1:.2f}  {2:.2f}'.format(d[6], d[7], d[8]))
    xp.xprint('Output:')
    xp.xprint('   GSM                x             y             z           MLT')
    xp.xprint('   SpacePy:    {0:12.8f}  {1:12.8f}  {2:12.8f}  {3:11.8f}' \
        .format(spacepy[0], spacepy[1], spacepy[2], spacepy[3]))
    xp.xprint('   Geopack:    {0:12.8f}  {1:12.8f}  {2:12.8f}  {3:11.8f}' \
        .format(geopack[0], geopack[1], geopack[2], geopack[3]))
    xp.xprint('   cxform:     {0:12.8f}  {1:12.8f}  {2:12.8f}  {3:11.8f}' \
        .format(cxform[0], cxform[1], cxform[2], cxform[3]))
    xp.xprint('   SSCWeb:     {0:12.8f}  {1:12.8f}  {2:12.8f}  {3:11.8f}' \
        .format(sscweb[k][0], sscweb[k][1], sscweb[k][2], sscweb[k][3]))
    #print('   Diff:       {0:.3f} {1:.3f} {2:.3f} {3:.3f}' \
    #    .format(sscweb[k][0] - spacepy[0],sscweb[k][1] - spacepy[1],
    #            sscweb[k][2] - spacepy[2],sscweb[k][3] - spacepy[3]))
    k = k + 1
