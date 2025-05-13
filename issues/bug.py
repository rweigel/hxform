import datetime
import numpy as np

from geopack import geopack
from hxform import hxform as hx

# TODO: What was issue?
time = (2003, 11, 19, 21, 36, 0)

ut = (datetime.datetime(*time)-datetime.datetime(1970,1,1)).total_seconds()
ps = geopack.recalc(ut)
print("Dipole tilt: {}".format(ps*180.0/np.pi))
#print(geopack.dip(-5.1,0.3,2.8))

dipole = hx.transform(np.array([0., 0., 1.]), time, 'MAG', 'GSM', 'car', 'sph')
print("dipole = {}".format(dipole))
print("Dipole tilt: {}".format(dipole[1]-90))

subsol_pt = hx.transform(np.array([0., 0., 1.]), time, 'MAG', 'GSM')
print("subsolar point = {}".format(subsol_pt))
phi = np.arctan2(subsol_pt[1], subsol_pt[0])
print("Subsolar angle = {}".format(phi*180/np.pi))