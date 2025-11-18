import astropy.coordinates
import sunpy.coordinates
# sunpy.coordinates is not used directly, but needed to register the frames.

import numpy
def transform(x, y, z, frame_in, frame_out, one):
  kwargs = {
    "x": x*one,
    "y": y*one,
    "z": z*one,
    "frame": frame_in,
    "obstime": '2010-12-30T00:00:00',
    "representation_type": 'cartesian'
  }
  coord = astropy.coordinates.SkyCoord(**kwargs)
  coord = coord.transform_to(frame_out).cartesian/one
  return coord.xyz.decompose().value.transpose()

# See https://github.com/sunpy/sunpy/pull/7530 for motivation for using
# different values for one. Seems not possible to transform dimensionless
# vectors in SunPy.
frames = [
  "geocentricearthequatorial",
  "geocentricsolarecliptic",
  "geocentricsolarmagnetospheric",
  "itrs",
  "solarmagnetic",
  "geomagnetic"
]

print(sunpy.__version__)
for frame_in in frames:
  for frame_out in frames:
    print(40*"-")
    for one in [astropy.constants.R_earth, astropy.units.m]:
      print(f"frame_in = {frame_in} | frame_out = {frame_out} | one = {one}")

      v = numpy.array([1.0, 1.0, 1.0])
      vt = transform(v[0], v[1], v[2], frame_in, frame_out, one)

      # Build transform matrix
      c1 = transform(1.0, 0.0, 0.0, frame_in, frame_out, one)
      c2 = transform(0.0, 1.0, 0.0, frame_in, frame_out, one)
      c3 = transform(0.0, 0.0, 1.0, frame_in, frame_out, one)
      matrix = numpy.column_stack([c1, c2, c3])

      diff = vt - numpy.dot(matrix, v)
      if numpy.all(numpy.abs(diff) < 1e-15):
        # cxform, Geopack-08, SpacePy, PySPEDAS, SpiceyPy, SSCWeb all pass.
        print("Pass")
      else:
        print("Fail")
        print(diff)

