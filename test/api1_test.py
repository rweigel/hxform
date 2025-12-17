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
ts = f"{t[0]:04d}-{t[1]:02d}-{t[2]:02d}T{t[3]:02d}:{t[4]:02d}:{t[5]:02d}"

match = "v must be a list, tuple, or np.ndarray"
assert_raises(ValueError, hxform.transform, [None, t, 'GEO', 'GSE'], match=match)

match = "v must be a list, tuple, or np.ndarray"
assert_raises(ValueError, hxform.transform, [{}, t, 'GEO', 'GSE'], match=match)

match = "v must be a list, tuple, or np.ndarray"
assert_raises(ValueError, hxform.transform, [{"a": 1}, t, 'GEO', 'GSE'], match=match)

match = "v is an empty list"
assert_raises(ValueError, hxform.transform, [[], t, 'GEO', 'GSE'], match=match)

match = "v is an empty tuple"
assert_raises(ValueError, hxform.transform, [(), t, 'GEO', 'GSE'], match=match)

match = "v may not be a scalar (0-dim ndarray)"
assert_raises(ValueError, hxform.transform, [np.array(1), t, 'GEO', 'GSE'], match=match)

match = "utilrsw.np.components2matrix(v) raised"
assert_raises(ValueError, hxform.transform, [v[0:2], t, 'GEO', 'GSE'], match=match)

match = "utilrsw.np.components2matrix(v) raised"
assert_raises(ValueError, hxform.transform, [np.array(v[0:2]), t, 'GEO', 'GSE'], match=match)

match = "utilrsw.np.components2matrix(v) raised"
vx = [1, 'a', 3]
assert_raises(ValueError, hxform.transform, [vx, t, 'GEO', 'GSE'], match=match)

match = "utilrsw.np.components2matrix(v) raised"
vx = ['a', 1, 2]
assert_raises(ValueError, hxform.transform, [vx, t, 'GEO', 'GSE'], match=match)

match = "utilrsw.np.components2matrix(v) raised"
vx = [{}, 1, 3]
assert_raises(ValueError, hxform.transform, [vx, t, 'GEO', 'GSE'], match=match)

tsx = '2009-02-02 11:01:01'
match = "datetime.datetime.strptime('2009-02-02 11:01:01', '%Y-%m-%dT%H:%M:%S') failed"
assert_raises(ValueError, hxform.transform, [v, tsx, 'GEO', 'GSE'], match=match)

match = "datetime.datetime.strptime('2009', '%Y-%m-%dT%H:%M:%S') failed"
assert_raises(ValueError, hxform.transform, [None, ts[0:4], 'GEO', 'GSE'], match=match)

match = "If number of vectors and number of times are both > 1, they must be equal"
assert_raises(ValueError, hxform.transform, [[v, v], [t, t, t], 'GEO', 'GSE'], match=match)

match = "If number of vectors and number of times are both > 1, they must be equal"
assert_raises(ValueError, hxform.transform, [[v, v], [ts, ts, ts], 'GEO', 'GSE'], match=match)

match = "If time is a list, tuple, or ndarray of strings, all elements must be strings"
assert_raises(ValueError, hxform.transform, [[v, v], [ts, t], 'GEO', 'GSE'], match=match)

match = "When represented as a list or tuple of ints, each time must have at least 3 elements"
assert_raises(ValueError, hxform.transform, [v, t[0:2], 'GEO', 'GSE'], match=match)

match = "When represented as a list or tuple of ints, each time must be given by 3 to 6 ints"
assert_raises(ValueError, hxform.transform, [v, t + [0], 'GEO', 'GSE'], match=match)

match = "Invalid time input"
assert_raises(ValueError, hxform.transform, [v, [[t]], 'GEO', 'GSE'], match=match)
assert_raises(ValueError, hxform.transform, [v, [[ts]], 'GEO', 'GSE'], match=match)
assert_raises(ValueError, hxform.transform, [v, [None], 'GEO', 'GSE'], match=match)

match = "time must be a str, list, tuple, or np.ndarray"
assert_raises(ValueError, hxform.transform, [v, None, 'GEO', 'GSE'], match=match)

libs = hxform.libs()
for lib in libs:

  coda = ""
  if skip_sscweb and lib == 'sscweb':
    coda = "(skipping b/c skip_sscweb=True)"
  hxform.xprint(f"Testing {lib} {coda}")
  if skip_sscweb and lib == 'sscweb':
    continue

  # Single vector, single timestamp
  vt = trans(v, t, lib)
  assert(isinstance(vt, list))  # Outer type of v is preserved
  assert(len(vt) == 3)
  assert(all_float(vt))         # v is ints, but output is float

  vt_s = trans(v, ts, lib)
  assert(isinstance(vt_s, list))  # Outer type of v is preserved
  assert(len(vt_s) == 3)
  assert(all_float(vt))           # v is ints, but vt is float
  assert(vt == vt_s)              # Outputs match when time formats differ


  vt = trans(np.array(v), np.array(t), lib)
  assert(isinstance(vt, np.ndarray))    # Outer type of v is preserved
  assert(vt.dtype == np.double)         # v is int64, but vt is np.double
  assert(vt.shape == (3, ))             # Shape matches input

  vt_s = trans(np.array(v), np.array(ts), lib)
  assert(isinstance(vt, np.ndarray))    # Outer type of input is preserved
  assert(vt.dtype == np.double)         # Input is int64, but output is np.double
  assert(vt.shape == (3, ))             # Shape matches input
  assert(np.all(vt == vt_s))            # Outputs match when time formats differ


  vt = trans(np.array(v), t, lib)
  assert(isinstance(vt, np.ndarray)) # Outer type of input is preserved
  assert(vt.dtype == np.double)      # Input is ints, but output is np.double
  assert(vt.shape == (3, ))

  vt_s = trans(np.array(v), ts, lib)
  assert(isinstance(vt_s, np.ndarray)) # Outer type of input is preserved
  assert(vt_s.dtype == np.double)      # Input is ints, but output is np.double
  assert(vt_s.shape == (3, ))
  assert(np.all(vt == vt_s))           # Outputs match when time formats differ


  vt = trans(v, np.array(t), lib)
  assert(isinstance(vt, list))       # Outer type of input is preserved
  assert(len(vt) == 3)
  assert(all_float(vt))              # v is ints, but vt is float

  vt_s = trans(v, np.array(ts), lib)
  assert(isinstance(vt, list))        # Outer type of input is preserved
  assert(len(vt_s) == 3)
  assert(len(vt_s) == 3)
  assert(all_float(vt_s))             # v is ints, but vt is float
  assert(vt == vt_s)                  # Outputs match when time formats differ


  vt = trans([v], [t], lib)
  assert(isinstance(vt, list))      # Outer type of input is preserved
  assert(isinstance(vt[0], list))   # Inner type of input is preserved
  assert(len(vt) == 1)
  assert(len(vt[0]) == 3)
  assert(all_float(vt))

  vt_s = trans([v], [ts], lib)
  assert(isinstance(vt, list))      # Outer type of input is preserved
  assert(isinstance(vt[0], list))   # Inner type of input is preserved
  assert(len(vt_s) == 1)
  assert(len(vt_s[0]) == 3)
  assert(all_float(vt_s))
  assert(vt == vt_s)                # Outputs match when time formats differ


  vt = trans([np.array(v)], t, lib)
  assert(isinstance(vt, list))           # Outer type of input is preserved
  assert(len(vt) == 1)
  assert(isinstance(vt[0], np.ndarray))  # Inner type of input is preserved
  assert(vt[0].shape == (3, ))

  vt_s = trans([np.array(v)], ts, lib)
  assert(isinstance(vt_s, list))           # Outer type of input is preserved
  assert(len(vt_s) == 1)
  assert(isinstance(vt_s[0], np.ndarray))  # Inner type of input is preserved
  assert(vt_s[0].shape == (3, ))
  assert(np.all(vt[0] == vt_s[0]))         # Outputs match when time formats differ


  vt = trans([v], np.array(t), lib)
  assert(isinstance(vt, list))           # Outer type of input is preserved
  assert(len(vt) == 1)
  assert(isinstance(vt[0], list))        # Inner type of input is preserved
  assert(len(vt[0]) == 3)

  vt_s = trans([v], np.array(ts), lib)
  assert(isinstance(vt_s, list))           # Outer type of input is preserved
  assert(len(vt_s) == 1)
  assert(isinstance(vt_s[0], list))        # Inner type of input is preserved
  assert(len(vt_s[0]) == 3)
  assert(vt == vt_s)


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

  vt1 = trans(v, [t, t], lib)
  vt2 = trans(v, [ts, ts], lib)
  assert(vt1 == vt2)

  vt1 = trans(v, np.array([ts, ts]), lib)
  vt2 = trans(v, [t, t], lib)
  assert(vt1 == vt2)

