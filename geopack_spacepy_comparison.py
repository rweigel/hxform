"""
hxtransform has not been compiled, execute
$ f2py -c Geopack-2008_dp.for T96_01.for -m hxtransform

This script produces 2 files
    1. error analysis file
    2. time analysis

you may vary arr_size to change the size of the arrays.


"""
import os
import sys
import numpy as np
import hxtransform as hx
import cxtransform as cx
from datetime import datetime
import time as tm

# Tsyganenko uses Geocentric-Solar Wind GSW instead of GSM. GSW has the positive
# x-axis pointing antiparallel to the solar wind. He chose this coordinate system
# in order to simplify his code. GSW becomes identical to GSM when the solar Wind
# velocity travels only in the negative x-axis direction, Tsyganenko chose
#  <-400,0,0> SMGSW_08
debug = False

arr_size = 10

X1 = np.arange(0,arr_size,dtype=np.float64)
Y1 = np.arange(0,arr_size,dtype=np.float64)
Z1 = np.arange(0,arr_size,dtype=np.float64)

trans = []
csys = ['GSM','SM']#,'GSE','MAG','GEO','GEI']
for c1 in csys:
    for c2 in csys:
        if c1 != c2:
            trans.append('to'.join((c1,c2)))

ferr = open('error_analysis.txt','w')
ftime = open('time_analysis.txt','w')

# allowable timeframe 1965-2020
# default is np.int64  but neet it to be np.int32 to match fortran INTEGER*4
# time is [YYYY,DDD,HH,MM,SS]
dtime = np.array((1997,1,0,0,0),dtype=np.int32)
for t in trans:
    print('\n    new\n   ')
    geo_start = tm.time()
    X2, Y2, Z2 = hx.transform(X1,Y1,Z1,t,dtime)
    # print(X2)
    geo_end = tm.time()

    time = "{0:04} {1:03d} {2:02d} {3:02d} {4:02d}".format(dtime[0],
                                                           dtime[1], dtime[2],
                                                           dtime[3], dtime[4])
    time = datetime.strptime(time, "%Y %j %H %M %S")
    time = [time.year,time.month, time.day, time.hour,time.minute,time.second]

    initial, final = t.split('to')

    sp_start = tm.time()
    spcoord = cx.transform(np.column_stack([X1,Y1,Z1]), time, initial, final)
    # print(spcoord[:,0])
    # print(X1,Y1,Z1)
    sp_end = tm.time()

    geocoord = np.column_stack([X2,Y2,Z2])
    magdiff = np.linalg.norm(geocoord-spcoord)/np.linalg.norm(spcoord)

    geo_t_diff = geo_end - geo_start
    sp_t_diff = sp_end - sp_start

    ferr.write("{} : {}\n".format(t, magdiff))
    ftime.write("{}: \n".format(t ))
    ftime.write("    spacepy: {} seconds\n".format(sp_t_diff))
    ftime.write("    geopack: {} seconds\n".format(geo_t_diff))
    ftime.write("    %difference: {} % \n".format( abs(sp_t_diff-geo_t_diff)*100 / ((sp_t_diff+geo_t_diff)/2)))
    # break
ftime.close()
ferr.close()
