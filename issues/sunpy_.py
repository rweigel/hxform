#import sunpy
from sunpy.coordinates import get_horizons_coord
#import sunpy
#sunpy.log.setLevel('DEBUG')
sc_id = -92 # ACE

data_jpl = get_horizons_coord(sc_id, "2021-11-25T00:00:00").geocentricsolarecliptic.cartesian.xyz.to('km')
print(data_jpl)
print(data_jpl.value/6378.16)

sc_id = -78 # DSCOVR
data_jpl = get_horizons_coord(sc_id, "2021-11-25T00:00:00").geocentricsolarecliptic.cartesian.xyz.to('km')
print(data_jpl)
print(data_jpl.value/6378.16)
