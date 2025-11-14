'''
API Test 3
Test that conversions between cartesian and spherical by hxform are consistent
with when doing the transform outside of hxform.
'''

import numpy
import hxform

skip_sscweb = True
skip_pyspedas = True

time  = [2010, 12, 30, 0, 0, 0]
input = numpy.array([1., 1., 1.])

libs = hxform.libs()

for lib in libs:
  if skip_sscweb and lib == 'sscweb':
    continue

  #hxform.xprint(f'Testing {lib}')


  output1 = hxform.transform(input, time, 'GEO', 'GSE', 'car', 'car', lib=lib)
  output1 = hxform.car2sph(output1)
  output2 = hxform.transform(input, time, 'GEO', 'GSE', 'car', 'sph', lib=lib)

  assert(numpy.all(output1 == output2))


  output1 = hxform.transform(input, time, 'GEO', 'GSE', 'car', 'sph', lib=lib)
  output2 = hxform.transform(input, time, 'GEO', 'GSE', 'car', 'car', lib=lib)
  output2 = hxform.car2sph(output2)

  assert(numpy.all(output1 == output2))


  output1 = hxform.transform(input, time, 'GEO', 'GSE', 'car', 'car', lib=lib)
  input2 = hxform.car2sph(input)
  output2 = hxform.transform(input2, time, 'GEO', 'GSE', 'sph', 'car', lib=lib)

  assert(numpy.all(numpy.abs(output1 - output2) < 1e-15))


  input1 = hxform.car2sph(input)
  output1 = hxform.transform(input1, time, 'GEO', 'GSE', 'sph', 'sph', lib=lib)
  output1 = hxform.sph2car(output1)

  output2 = hxform.transform(input1, time, 'GEO', 'GSE', 'sph', 'car', lib=lib)

  assert(numpy.all(numpy.abs(output1 - output2) < 1e-15))
