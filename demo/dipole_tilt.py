import datetime
import numpy as np
import hxform

from geopack import geopack

# TODO: What was issue this was written for?
time = (2003, 11, 19, 21, 36, 0)
input = np.array([0., 0., 1.])
hxform.xprint(f"Time:  {time}")
hxform.xprint(f"Input: {input}")

ut = (datetime.datetime(*time)-datetime.datetime(1970,1,1)).total_seconds()
ps = geopack.recalc(ut)
hxform.xprint("Geopack dipole tilt: {}".format(ps*180.0/np.pi))

output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'car', lib='spacepy')
hxform.xprint(f"output (cartesian) = {output}")
hxform.xprint(f"output (spherical) = {hxform.car2sph(*output)}")
output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'sph', lib='spacepy')
hxform.xprint(f"output (spherical) = {output}")
hxform.xprint(f"dipole tilt = {output[1]-90}")

#hxform.xprint("Dipole tilt: {}".format(dipole[1]-90))

output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'car', lib='geopack_08_dp')
hxform.xprint(f"output (cartesian) = {output}")
hxform.xprint(f"output (spherical) = {hxform.car2sph(*output)}")
output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'sph', lib='geopack_08_dp')
hxform.xprint(f"output (spherical) = {output}")
hxform.xprint(f"dipole tilt = {output[1]-90}")

#subsol_pt = hxform.transform(input, time, 'MAG', 'GSM')
#phi = np.arctan2(subsol_pt[1], subsol_pt[0])
#hxform.xprint("subsolar point = {}".format(subsol_pt))
#hxform.xprint("Subsolar angle = {}".format(phi*180/np.pi))