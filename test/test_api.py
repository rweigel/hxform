# API Demo
import numpy as np
from hxform import hxform as hx
from hxform import xprint

time1  = [2010, 12, 30, 0, 0, 0]
input1 = [1., 1., 1.]
time2  = [time1, time1]
input2 = [input1, input1]

# With
# csys_in = 'GSM'
# csys_out = 'GSE'
# geopack_08_dp gives error for
#  assert(output12[0] == output11)
#  assert(output12[1] == output11)


libs = hx.known_libs(info=False)

for lib in libs:
  if lib == 'geopack_08_dp':
    continue
  if lib != 'sscweb':
    continue

  for t1 in hx.known_transforms(lib):
    for t2 in hx.known_transforms(lib):

      csys_in = t1
      csys_out = t2

      xprint(f'{lib} {t1} to {t2}')

      # Single time, single vector
      output11 = hx.transform(input1, time1, csys_in, csys_out, lib=lib)
      print(output11)
      assert(output11 == output11)

      # Multiple times, single vector
      output12 = hx.transform(input1, time2, csys_in, csys_out, lib=lib)
      print(output12)
      assert(output12[0] == output11)
      assert(output12[1] == output11)

      # Single time, multiple vectors
      output21 = hx.transform(input2, time1, csys_in, csys_out, lib=lib)
      print(output21)
      assert(output21[0] == output11)
      assert(output21[1] == output11)

      # Multiple times, multiple vectors
      output22 = hx.transform(input2, time2, csys_in, csys_out, lib=lib)
      print(output22)
      assert(output22[0] == output11)
      assert(output22[1] == output11)
