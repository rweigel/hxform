import datetime

import numpy

import hxform

from geopack import geopack

#time = (2003, 11, 19, 21, 36, 0)
time = (2010, 12, 21, 0, 0, 0)
#time = (1970, 1, 1, 0, 0, 0)
input = [0., 0., 1.]
hxform.xprint(f"Time:  {time}")

ut = (datetime.datetime(*time)-datetime.datetime(1970,1,1)).total_seconds()
ps = geopack.recalc(ut)
hxform.xprint("Geopack dipole tilt using recalc():  {}".format(ps*180.0/numpy.pi))

# Why does this differ from the recalc() value 0.0001 degrees? For double
# precision, it should be same to ~1e-15 degrees.
#output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'car', lib='geopack_08_dp')
#hxform.xprint(f"output (cartesian) = {output}")
#hxform.xprint(f"output (spherical) = {hxform.car2sph(*output)}")
output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'sph', lib='geopack_08_dp')
#hxform.xprint(f"output (spherical) = {output}")
hxform.xprint(f"Geopack dipole tilt using transform: {output[1]-90}")

#output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'car', lib='spacepy')
#hxform.xprint(f"output (cartesian) = {output}")
#hxform.xprint(f"output (spherical) = {hxform.car2sph(*output)}")
output = hxform.transform(input, time, 'MAG', 'GSM', 'car', 'sph', lib='spacepy')
#hxform.xprint(f"output (spherical) = {output}")
hxform.xprint(f"SpacePy dipole tilt using transform: {output[1]-90}")

