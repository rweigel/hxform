import os
import spiceypy as sp;

from hxform import xprint # print to console and spiceypy_demo.log

kernels = ['naif0012.tls', 'rbsp_general011.tf', 'de440s.bsp', 'pck00011.tpc']
for kernel in kernels:
  rel_path = os.path.join('spiceypy', 'kernels', kernel)
  kernel_file = os.path.join(os.path.dirname(__file__), rel_path)
  xprint(f"Furnishing: {rel_path}")
  sp.furnsh(kernel_file)

# Ephemeris time (TDB seconds past J2000.0) is the time currency 
# of SPICE.  Start by converting a UTC epoch of interest into it.
et = sp.str2et("2013-AUG-10 12:00:00.000")

# Now compute the rotation matrix from GEO to GSE at this time.
geo2gse = sp.pxform("GEO", "GSE", et)

# And the rotation matrix from GSE to GSM.
gse2gsm = sp.pxform("GSE", "GSM", et)

# Apply the rotations to the vector of interest.
vec_geo = [1,0,0]
xprint(f' vec_geo: {vec_geo}')
vec_gse = sp.mxv(geo2gse, vec_geo)
xprint(f' vec_gse: {vec_gse}')
vec_gsm = sp.mxv(gse2gsm, vec_gse)
xprint(f' vec_gsm: {vec_gsm}')
# Unload SPICE kernels
sp.kclear()