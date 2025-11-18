'''
API Test 1
Test that output type matches input type.
Test that output shape is correct.
'''

import numpy as np
import hxform

skip_sscweb = True

def trans(v, t, lib):
  # lib should not matter.
  return hxform.transform(v, t, 'GEO', 'GSE', lib=lib)

def all_float(x):
  for i in range(len(x)):
    if isinstance(x[i], list):
      if not all_float(x[i]):
        return False
    else:
      if not isinstance(x[i], float):
        return False
  return True


v = [1, 2, 3]
t = [2009, 2, 2, 11, 1, 1]

libs = hxform.libs()
for lib in libs:
  coda = ""
  if skip_sscweb and lib == 'sscweb':
    coda = "(skipping b/c skip_sscweb=True)"
  hxform.xprint(f"Testing {lib} {coda}")
  if skip_sscweb and lib == 'sscweb':
    continue

  vt = trans(v, t, lib)
  assert(isinstance(vt, list))  # Outer type of v is preserved
  assert(len(vt) == 3)
  assert(all_float(vt)) # Input is ints, but output is float

  vt = trans(np.array(v), np.array(t), lib)
  assert(isinstance(vt, np.ndarray))    # Outer type of v is preserved
  assert(vt.dtype == np.double)         # Input is ints, but output is np.double
  assert(vt.shape == (3, ))

  vt = trans(np.array(v), t, lib)
  assert(isinstance(vt, np.ndarray)) # Outer type of v is preserved
  assert(vt.dtype == np.double)
  assert(vt.shape == (3, ))

  vt = trans(v, np.array(t), lib)
  assert(isinstance(vt, list)) # Outer type of v is preserved
  assert(len(vt) == 3)

  vt = trans([v],[t], lib)
  assert(isinstance(vt, list))
  assert(isinstance(vt[0], list))
  assert(len(vt) == 1)
  assert(len(vt[0]) == 3)
  assert(all_float(vt))

  vp = trans([np.array(v)], t, lib)
  assert(isinstance(vp, list)) # Outer type of v is preserved
  assert(len(vp) == 1)
  assert(isinstance(vp[0], np.ndarray))
  assert(vp[0].shape == (3, ))

  vp = trans([v], np.array(t), lib)
  assert(isinstance(vp, list)) # Outer type of v is preserved
  assert(len(vp) == 1)
  assert(isinstance(vp[0], list))
  assert(len(vp[0]) == 3)

  vt = trans(np.array([v]), np.array([t]), lib)
  assert(isinstance(vt, np.ndarray))
  assert(vt.shape == (1, 3))
  assert(isinstance(vt[0, 0], np.double))

  vt = trans([v, v],[t, t], lib)
  assert(isinstance(vt, list))
  assert(isinstance(vt[0], list))
  assert(isinstance(vt[1], list))
  assert(len(vt) == 2)
  assert(len(vt[0]) == 3)
  assert(len(vt[1]) == 3)
  assert(all_float(vt))

  vt = trans(np.array([v, v]), np.array([t,t]), lib)
  assert(isinstance(vt, np.ndarray))
  assert(isinstance(vt[0], np.ndarray))
  assert(isinstance(vt[1], np.ndarray))
  assert(vt.shape == (2, 3))
  assert(vt[0].dtype == np.double)
  assert(vt[1].dtype == np.double)

  vt = trans(np.array(v), np.array([t,t]), lib)
  assert(vt.shape == (2, 3))

  vt = trans(np.array([v, v]), t, lib)
  assert(vt.shape == (2, 3))

  vt = trans(v, [t, t], lib)
  assert(isinstance(vt, list))
  assert(len(vt) == 2)
  assert(isinstance(vt[0], list))
  assert(isinstance(vt[1], list))

  vt = trans([v, v], t, lib)
  assert(isinstance(vt, list))
  assert(len(vt) == 2)
  assert(isinstance(vt[0], list))
  assert(isinstance(vt[1], list))
  assert(len(vt[0]) == 3)
  assert(len(vt[1]) == 3)
