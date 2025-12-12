'''
API Test 1
Test that output outer and inner type matches that in input.
Test that output lengths/shapes match input lengths/shapes.
'''

import numpy as np
import hxform

from utilrsw.test.assert_raises import assert_raises

skip_sscweb = True

def trans(v, t, lib):
  # lib, frame_in, and frame_out should not matter.
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

  match = "utilrsw.np.components2matrix(v) raised"
  assert_raises(ValueError, hxform.transform, [None, t, 'GEO', 'GSE'], match=match)

  match = "utilrsw.np.components2matrix(v) raised"
  assert_raises(ValueError, hxform.transform, [v[0:2], t, 'GEO', 'GSE'], match=match)

  match = "utilrsw.np.components2matrix(v) raised"
  assert_raises(ValueError, hxform.transform, [{}, t, 'GEO', 'GSE'], match=match)

  # Single vector, single timestamp
  vt = trans(v, t, lib)
  assert(isinstance(vt, list))  # Outer type of input is preserved
  assert(len(vt) == 3)
  assert(all_float(vt))         # Input is ints, but output is float

  vt = trans(np.array(v), np.array(t), lib)
  assert(isinstance(vt, np.ndarray))    # Outer type of input is preserved
  assert(vt.dtype == np.double)         # Input is ints, but output is np.double
  assert(vt.shape == (3, ))

  vt = trans(np.array(v), t, lib)
  assert(isinstance(vt, np.ndarray)) # Outer type of input is preserved
  assert(vt.dtype == np.double)      # Input is ints, but output is np.double
  assert(vt.shape == (3, ))

  vt = trans(v, np.array(t), lib)
  assert(isinstance(vt, list))       # Outer type of input is preserved
  assert(len(vt) == 3)

  vt = trans([v], [t], lib)
  assert(isinstance(vt, list))      # Outer type of input is preserved
  assert(isinstance(vt[0], list))   # Inner type of input is preserved
  assert(len(vt) == 1)
  assert(len(vt[0]) == 3)
  assert(all_float(vt))

  vt = trans([np.array(v)], t, lib)
  assert(isinstance(vt, list))           # Outer type of input is preserved
  assert(len(vt) == 1)
  assert(isinstance(vt[0], np.ndarray))  # Inner type of input is preserved
  assert(vt[0].shape == (3, ))

  vt = trans([v], np.array(t), lib)
  assert(isinstance(vt, list))           # Outer type of input is preserved
  assert(len(vt) == 1)
  assert(isinstance(vt[0], list))        # Inner type of input is preserved
  assert(len(vt[0]) == 3)

  vt = trans(np.array([v]), np.array([t]), lib)
  assert(isinstance(vt, np.ndarray))      # Outer type of input is preserved
  assert(vt.shape == (1, 3))
  assert(isinstance(vt[0, 0], np.double)) # Input is ints, but output is np.double


  # Single vector, multiple timestamps
  vt = trans(v, [t, t], lib)
  assert(isinstance(vt, list))    # Outer type of input is preserved
  assert(len(vt) == 2)            # Number of output items matches number of timestamps
  assert(isinstance(vt[0], list)) # Inner type of input is preserved
  assert(isinstance(vt[1], list)) # Inner type of input is preserved

  vt = trans(np.array(v), [t, t], lib)
  assert(isinstance(vt, np.ndarray)) # Outer type of input is preserved
  assert(vt.shape == (2, 3))         # Number of rows matches number of timestamps

  vt = trans(v, np.array([t, t]), lib)
  assert(isinstance(vt, list))    # Outer type of input is preserved
  assert(len(vt) == 2)            # Number of output items matches number of timestamps
  assert(isinstance(vt[0], list)) # Inner type of input is preserved
  assert(isinstance(vt[1], list)) # Inner type of input is preserved

  vt = trans(np.array(v), np.array([t, t]), lib)
  assert(isinstance(vt, np.ndarray))
  assert(vt.shape == (2, 3))      # Number of rows matches number of timestamps


  # Multiple vectors, single timestamp
  vt = trans([v, v], t, lib)
  assert(isinstance(vt, list))    # Outer type of input is preserved
  assert(len(vt) == 2)
  assert(isinstance(vt[0], list)) # Inner type of input is preserved
  assert(isinstance(vt[1], list)) # Inner type of input is preserved
  assert(len(vt[0]) == 3)
  assert(len(vt[1]) == 3)

  vt = trans([v, v], np.array(t), lib)
  assert(isinstance(vt, list))    # Outer type of input is preserved
  assert(len(vt) == 2)
  assert(isinstance(vt[0], list)) # Inner type of input is preserved
  assert(isinstance(vt[1], list)) # Inner type of input is preserved
  assert(len(vt[0]) == 3)
  assert(len(vt[1]) == 3)

  vt = trans(np.array([v, v]), t, lib)
  assert(isinstance(vt, np.ndarray))
  assert(vt.shape == (2, 3))  # Number of rows matches number of input rows

  vt = trans(np.array([v, v]), np.array(t), lib)
  assert(isinstance(vt, np.ndarray))
  assert(vt.shape == (2, 3))  # Number of rows matches number of input rows


  # Multiple vectors, multiple timestamps
  vt = trans([v, v], [t, t], lib)
  assert(isinstance(vt, list))    # Outer type of input is preserved
  assert(isinstance(vt[0], list)) # Inner type of input is preserved
  assert(isinstance(vt[1], list)) # Inner type of input is preserved
  assert(len(vt) == 2)
  assert(len(vt[0]) == 3)
  assert(len(vt[1]) == 3)
  assert(all_float(vt))           # Input is ints, but output is float

  vt = trans([v, v], np.array([t, t]), lib)
  assert(isinstance(vt, list))    # Outer type of input is preserved
  assert(isinstance(vt[0], list)) # Inner type of input is preserved
  assert(isinstance(vt[1], list)) # Inner type of input is preserved
  assert(len(vt) == 2)
  assert(len(vt[0]) == 3)
  assert(len(vt[1]) == 3)
  assert(all_float(vt))           # Input is ints, but output is float

  vt = trans(np.array([v, v]), [t, t], lib)
  assert(isinstance(vt, np.ndarray))
  assert(vt.shape == (2, 3))  # Number of rows matches number of input rows

  vt = trans(np.array([v, v]), np.array([t, t]), lib)
  assert(isinstance(vt, np.ndarray))    # Outer type of input is preserved
  assert(vt.shape == (2, 3))            # Number of rows matches number of timestamps
  assert(vt[0].dtype == np.double)
  assert(vt[1].dtype == np.double)
