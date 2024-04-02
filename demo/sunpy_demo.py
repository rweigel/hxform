# This adds the coordinates and units method to top-level astropy package.
import astropy.coordinates

# Add SunPy coordinate frames in astropy.coordinates. Editor may note that
# import is not used, but it is.
import sunpy.coordinates

from hxform import xprint # print to console and spiceypy_demo.log

frames = astropy.coordinates.frame_transform_graph.get_names()
frames.sort()
xprint("Available frames:")
for frame in frames:
  xprint(" " + frame)

# https://datacenter.iers.org/eop.php
# https://datacenter.iers.org/data/9/finals2000A.all

name_map = {
  "GEI": "geocentricearthequatorial",
  "GSE": "geocentricsolarecliptic",
  "GSM": "geocentricsolarmagnetospheric",
  "GEO": "itrs",
  "SM":  "solarmagnetic",
  "MAG": "geomagnetic"
}

time = "2010-12-30T00:00"
initial = 'GSM'
final = 'GSE'
initial = 'GEI'
final = 'GSM'

R_E = astropy.constants.R_earth
# TODO: Use
#   R_E = astropy.units.m
# when https://github.com/sunpy/sunpy/pull/7530 is merged.
kwargs = {
  "x": [R_E,R_E],
  "y": [R_E,R_E],
  "z": [R_E,R_E],
  "frame": name_map[initial],
  "obstime": time,
  "representation_type": "cartesian"
}

coord_in = astropy.coordinates.SkyCoord(**kwargs)

coord_out = coord_in.transform_to(name_map[final]).cartesian/R_E

vout = coord_out.xyz.decompose()
xprint('Using sunpy directly')
xprint(vout.value.transpose())

import numpy as np
xtime = np.array([2010, 12, 30, 0, 0, 0])
v = np.array([[1,1,1],[1,1,1]])
from hxform import hxform as hx
vp = hx.transform(v, xtime, initial, final, lib='sunpy')
xprint('Using hxform wrapper')
xprint(vp)

