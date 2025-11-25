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

libs_only = []
#libs_only = ['sunpy'] # Uncomment to test only one or more specific libraries

raise_on_fail = False

def report(test_passed, diff, raise_on_fail, key, msg):
  if not test_passed:
    hxform.xprint(f"  FAIL for {msg}")
    hxform.xprint(f"  {diff}")
    if raise_on_fail:
      assert test_passed

    if key not in fails:
      fails[key] = []
    fails[key].append({'msg': msg, 'diff': diff})


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


      # Transform single vector at single timestamp
      msg = "single vector, single timestamp"
      try:
        vt = hxform.transform(v1, t1, f1, f2, lib=lib)
        matrix = hxform.matrix(t1, f1, f2, lib=lib)
        diff = vt - numpy.dot(matrix, v1)
        test_passed = numpy.all(numpy.abs(diff) < 1e-15)
        report(test_passed, diff, raise_on_fail, key, msg)
      except Exception as e:
        hxform.xprint(f"  FAIL on {msg}: {e}")
        if raise_on_fail:
          raise e

      # Transform two vectors at with same timestamp
      msg = "two vectors, single timestamp"
      try:
        vt = hxform.transform(v2, t1, f1, f2, lib=lib)
        matrix = hxform.matrix(t1, f1, f2, lib=lib)
        for i in range(v2.shape[0]):
          diff = vt[i, :] - numpy.dot(matrix, v2[i, :])
          report(test_passed, diff, raise_on_fail, key, f"{msg} vector {i}")
      except Exception as e:
        hxform.xprint(f"  FAIL on {msg}: {e}")
        if raise_on_fail:
          raise e

      # Transform two vectors at two timestamps
      msg = "two vectors, two timestamps"
      try:
        vt = hxform.transform(v2, t2, f1, f2, lib=lib)
        matrices = hxform.matrix(t2, f1, f2, lib=lib)
        for i in range(matrices.shape[0]):
          diff = vt[i, :] - numpy.dot(matrices[i, :, :], v2[i, :])
          report(test_passed, diff, raise_on_fail, key,  f"{msg} vector {i}")
      except Exception as e:
        hxform.xprint(f"  FAIL on {msg}: {e}")
        if raise_on_fail:
          raise e

if len(fails) != 0:
  import json
  fails = json.dumps(fails, indent=2, default=str)
  hxform.xprint(f"Some tests failed: \n{fails}")
  assert len(fails) == 0, "Some tests failed."
