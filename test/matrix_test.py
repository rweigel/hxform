"""
Test that transforming vector using dot(matrix(...), v) gives result consistent
with transform(v, ...).
"""
import numpy
import hxform

t1 = numpy.array([2010, 12, 30, 0, 0, 0])
v1 = numpy.array([1., 1., 1.])
v2 = numpy.array([v1, 2*v1])
t2 = numpy.array([t1, t1])

# SSCWeb returns results to only two decimal places, leading to large
# discrepancies when comparing matrix multiplication to transform results.
# TODO: Compute appropriate threshold for pass.
libs_exclude = ['sscweb']

libs_only = ['geopack_08_dp']
# If failure for SunPy, see https://github.com/sunpy/sunpy/issues/8406 and
# update setup.py to require version with fix.
#libs_only = ['sunpy'] # Uncomment to test only one or more specific libraries

raise_on_fail = False

def report(test_passed, diff, raise_on_fail, key, msg, error=None):
  if not test_passed:
    hxform.xprint(f"  {msg}")
    if error is not None:
      hxform.xprint(f"    FAIL. Error: {error}")
    else:
      error = "Difference exceeds threshold"
      hxform.xprint(f"    FAIL. diff = {diff}")

    if raise_on_fail:
      assert test_passed

    if key not in fails:
      fails[key] = []
    fails[key].append({'msg': msg, 'diff': diff, 'error': error})


fails = {}
libs = hxform.libs()
for lib in libs:

  if libs_only and lib not in libs_only:
    hxform.xprint(f"{lib}: Skipping tests because it is not in libs_only")
    continue

  if libs_exclude and lib in libs_exclude:
    hxform.xprint(f"{lib}: Skipping tests because it is in libs_exclude")
    continue

  for f1 in hxform.frames(lib):
    for f2 in hxform.frames(lib):

      key = f'{lib} {f1} to {f2}'
      hxform.xprint(key)

      if f1 == f2:
        continue

      if True:
        # Transform single vector at single timestamp
        msg = "matrix test: single vector, single timestamp"
        try:
          vt = hxform.transform(v1, t1, f1, f2, lib=lib)
          matrix = hxform.matrix(t1, f1, f2, lib=lib)
          diff = vt - numpy.dot(matrix, v1)
          test_passed = numpy.all(numpy.abs(diff) < 1e-15)
          report(test_passed, diff, raise_on_fail, key, msg)
        except Exception as e:
          print(e)
          report(False, None, raise_on_fail, key, msg, error=str(e))

        # Transform two vectors at with same timestamp
        msg = "matrix test: two vectors, single timestamp"
        try:
          vt = hxform.transform(v2, t1, f1, f2, lib=lib)
          matrix = hxform.matrix(t1, f1, f2, lib=lib)
          for i in range(v2.shape[0]):
            diff = vt[i, :] - numpy.dot(matrix, v2[i, :])
            test_passed = numpy.all(numpy.abs(diff) < 1e-15)
            report(test_passed, diff, raise_on_fail, key, f"{msg} vector {i}")
        except Exception as e:
          report(False, None, raise_on_fail, key, msg, error=str(e))



        # Transform two vectors at two timestamps
        msg = "matrix test: two vectors, two timestamps"
        try:
          vt = hxform.transform(v2, t2, f1, f2, lib=lib)
          matrices = hxform.matrix(t2, f1, f2, lib=lib)
          for i in range(matrices.shape[0]):
            diff = vt[i, :] - numpy.dot(matrices[i, :, :], v2[i, :])
            test_passed = numpy.all(numpy.abs(diff) < 1e-15)
            report(test_passed, diff, raise_on_fail, key,  f"{msg} vector {i}")
        except Exception as e:
          report(False, None, raise_on_fail, key, msg, error=str(e))


        msg = "car/sph test: test #1 on single vector, single timestamp"
        try:
          output1 = hxform.transform(v1, t1, f1, f2, 'car', lib=lib)
          output1 = hxform.car2sph(output1)
          output2 = hxform.transform(v1, t1, f1, f2, 'car', 'sph', lib=lib)
          diff = output1 - output2
          test_passed = numpy.all(numpy.abs(diff) < 1e-15)
          report(test_passed, diff, raise_on_fail, key, msg)
        except Exception as e:
          report(False, None, raise_on_fail, key, msg, error=str(e))


        msg = "car/sph test: test #2 on single vector, single timestamp"
        try:
          output1 = hxform.transform(v1, t1, f1, f2, 'car', 'sph', lib=lib)
          output2 = hxform.transform(v1, t1, f1, f2, 'car', 'car', lib=lib)
          output2 = hxform.car2sph(output2)
          diff = output1 - output2
          test_passed = numpy.all(numpy.abs(diff) < 1e-15)
          report(test_passed, diff, raise_on_fail, key, msg)
        except Exception as e:
          report(False, None, raise_on_fail, key, msg, error=str(e))


        msg = "car/sph test: test #3 on single vector, single timestamp"
        try:
          print(v1)
          print(t1)
          output1 = hxform.transform(v1, t1, f1, f2, 'car', 'car', lib=lib)
          print(v1)
          print(t1)
          input2 = hxform.car2sph(v1)
          print(v1)
          print(t1)
          output2 = hxform.transform(input2, t1, f1, f2, 'sph', 'car', lib=lib)
          diff = output1 - output2
          test_passed = numpy.all(numpy.abs(diff) < 1e-15)
          report(test_passed, diff, raise_on_fail, key, msg)
        except Exception as e:
          report(False, None, raise_on_fail, key, msg, error=str(e))
          exit()

      if False:

        msg = "car/sph test: test #4 on single vector, single timestamp"
        try:
          input1 = hxform.car2sph(v1)
          output1 = hxform.transform(input1, t1, f1, f2, 'sph', 'sph', lib=lib)
          output1 = hxform.sph2car(output1)
          output2 = hxform.transform(input1, t1, f1, f2, 'sph', 'car', lib=lib)
          diff = output1 - output2
          test_passed = numpy.all(numpy.abs(diff) < 1e-15)
          report(test_passed, diff, raise_on_fail, key, msg)
        except Exception as e:
          report(False, None, raise_on_fail, key, msg, error=str(e))


if len(fails) != 0:
  import json
  fails = json.dumps(fails, indent=2, default=str)
  hxform.xprint(f"Some tests failed: \n{fails}")
  assert len(fails) == 0, "Some tests failed."
