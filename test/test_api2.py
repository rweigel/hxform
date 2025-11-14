'''
API Test
For all transforms available in each library, test that the output is the
same for single and multiple times and vectors.
'''

import hxform as hx

skip_sscweb = True
skip_pyspedas = True

time1  = [2010, 12, 30, 0, 0, 0]
input1 = [1., 1., 1.]
time2  = [time1, time1]
input2 = [input1, input1]

libs = hx.libs()

for lib in libs:
  if skip_pyspedas and lib == 'pyspedas':
    print("pyspedas - Skipping due to https://github.com/spedas/pyspedas/issues/1274")
    continue

  if skip_sscweb and lib == 'sscweb':
    print("sscweb - Skipping b/c skip_sscweb=True")
    continue

  for t1 in hx.frames(lib):
    for t2 in hx.frames(lib):

      csys_in = t1
      csys_out = t2

      print(f'{lib} {t1} to {t2}')

      # Single time, single vector
      output11 = hx.transform(input1, time1, csys_in, csys_out, lib=lib)

      # Multiple times, single vector
      output12 = hx.transform(input1, time2, csys_in, csys_out, lib=lib)
      assert(output12[0] == output11)
      assert(output12[1] == output11)

      # Single time, multiple vectors
      output21 = hx.transform(input2, time1, csys_in, csys_out, lib=lib)
      assert(output21[0] == output11)
      assert(output21[1] == output11)

      # Multiple times, multiple vectors
      output22 = hx.transform(input2, time2, csys_in, csys_out, lib=lib)
      assert(output22[0] == output11)
      assert(output22[1] == output11)
