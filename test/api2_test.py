'''
API Test 2
For all transforms available in each library, test that the output is the
same for single and multiple times and vectors.
'''

import hxform

skip_sscweb = True
skip_pyspedas = True

time1  = [2010, 12, 30, 0, 0, 0]
input1 = [1., 1., 1.]
time2  = [time1, time1]
input2 = [input1, input1]

libs = hxform.libs()

for lib in libs:
  if skip_pyspedas and lib == 'pyspedas':
    print("pyspedas - Skipping due to https://github.com/spedas/pyspedas/issues/1274")
    continue

  if skip_sscweb and lib == 'sscweb':
    hxform.xprint("sscweb - Skipping b/c skip_sscweb=True")
    continue

  if lib != 'spacepy-irbem':
    continue

  for f1 in hxform.frames(lib):
    for f2 in hxform.frames(lib):

      frame_in = f1
      frame_out = f2

      hxform.xprint(f'{lib} {f1} to {f2}')

      # Single time, single vector
      #print(f"  Input: {input1} at {time1}")
      output11 = hxform.transform(input1, time1, frame_in, frame_out, lib=lib)

      # Multiple times, single vector
      #print(f"  Input: {input1} at {time2}")
      output12 = hxform.transform(input1, time2, frame_in, frame_out, lib=lib)
      assert(output12[0] == output11)
      assert(output12[1] == output11)

      # Single time, multiple vectors
      #print(f"  Input: {input2} at {time1}")
      output21 = hxform.transform(input2, time1, frame_in, frame_out, lib=lib)
      assert(output21[0] == output11)
      assert(output21[1] == output11)

      # Multiple times, multiple vectors
      #print(f"  Input: {input2} at {time2}")
      output22 = hxform.transform(input2, time2, frame_in, frame_out, lib=lib)
      assert(output22[0] == output11)
      assert(output22[1] == output11)
