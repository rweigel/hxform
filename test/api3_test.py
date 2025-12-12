import numpy

import hxform

from utilrsw.test.assert_raises import assert_raises

# lib, frame_in, and frame_out should not matter.

v = [1, 2, 3]
ti = [2009, 2, 2, 11, 1, 1]
ts = f"{ti[0]:04d}-{ti[1]:02d}-{ti[2]:02d}T{ti[3]:02d}:{ti[4]:02d}:{ti[5]:02d}"

vt1 = hxform.transform(v, ti, 'GEO', 'GSE')
vt2 = hxform.transform(v, ts, 'GEO', 'GSE')
assert(vt1 == vt2)

vt1 = hxform.transform(v, [ti], 'GEO', 'GSE')
vt2 = hxform.transform(v, [ts], 'GEO', 'GSE')
assert(vt1 == vt2)

vt1 = hxform.transform(v, [ti, ti], 'GEO', 'GSE')
vt2 = hxform.transform(v, [ts, ts], 'GEO', 'GSE')
assert(vt1 == vt2)

vt1 = hxform.transform(v, numpy.array([ts, ts]), 'GEO', 'GSE')
vt2 = hxform.transform(v, [ti, ti], 'GEO', 'GSE')
assert(vt1 == vt2)

tsx = '2009-02-02 11:01:01'
match = "datetime.datetime.strptime"
assert_raises(ValueError, hxform.transform, [v, tsx, 'GEO', 'GSE'], match=match)

match = "If number of vectors and number of times are both > 1, they must be equal"
assert_raises(ValueError, hxform.transform, [[v, v], [ti, ti, ti], 'GEO', 'GSE'], match=match)

match = "If number of vectors and number of times are both > 1, they must be equal"
assert_raises(ValueError, hxform.transform, [[v, v], [ts, ts, ts], 'GEO', 'GSE'], match=match)

match = "If time is a list, tuple, or ndarray of strings, all elements must be strings"
assert_raises(ValueError, hxform.transform, [[v, v], [ts, ti], 'GEO', 'GSE'], match=match)

match = "When represented as a list or tuple of ints, each time must have at least 3 elements"
assert_raises(ValueError, hxform.transform, [v, ti[0:2], 'GEO', 'GSE'], match=match)

match = "When represented as a list or tuple of ints, each time must be given by 3 to 6 ints"
assert_raises(ValueError, hxform.transform, [v, ti + [0], 'GEO', 'GSE'], match=match)

match = "Invalid time input"
assert_raises(ValueError, hxform.transform, [v, [[ti]], 'GEO', 'GSE'], match=match)
assert_raises(ValueError, hxform.transform, [v, [[ts]], 'GEO', 'GSE'], match=match)
assert_raises(ValueError, hxform.transform, [v, [None], 'GEO', 'GSE'], match=match)

match = "time must be a str, list, tuple, or np.ndarray"
assert_raises(ValueError, hxform.transform, [v, None, 'GEO', 'GSE'], match=match)
