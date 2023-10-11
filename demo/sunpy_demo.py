import astropy.units as u
from astropy.constants import R_earth
from astropy.coordinates import SkyCoord
import sunpy.coordinates  # registers GSE as a coordinate frame

# https://datacenter.iers.org/eop.php
# https://datacenter.iers.org/data/9/finals2000A.all
geo_coord = SkyCoord(x=R_earth, y=0*u.m, z=0*u.m,
                    frame='itrs', 
                    obstime='2013-08-10 12:00', representation_type='cartesian')

print(geo_coord.cartesian/R_earth)

print(geo_coord.transform_to('geocentricsolarecliptic').cartesian/R_earth)
