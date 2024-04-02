# API Test
# Test that output type matches input type.
# Test that output shape is correct.

import numpy as np
from hxform import hxform as hx

def trans(v,t):
  # lib should not matter.
  return hx.transform(v, t, 'GSM','MAG', lib='geopack_08_dp')


v = [1, 2, 3]
t = [2009, 2, 2, 11, 1, 1]

#############################################################################
vt = trans(v, t)
#print(v)
#print(vt)
assert(len(vt) == 3)
assert(type(vt) == type(v))  # Outer type is preserved
for i in range(3):
  assert(type(vt[i]) == float) # Input is ints, but output is always float type

vt = trans([v],[t])
#print([v])
#print(vt)
#print(type(vt))
assert(type(vt) == list)
assert(type(vt[0]) == list)
assert(len(vt) == 1)
assert(len(vt[0]) == 3)
for i in range(len(vt)):
  for j in range(3):
    assert(type(vt[i][j]) == float)

vt = trans([v,v],[t,t])
#print([v,v])
#print(vt)
assert(type(vt) == list)
assert(type(vt[0]) == list)
assert(type(vt[1]) == list)
assert(len(vt) == 2)
assert(len(vt[0]) == 3)
assert(len(vt[1]) == 3)
for i in range(len(vt)):
  for j in range(3):
    assert(type(vt[i][j]) == float)
#############################################################################


#############################################################################
vt = trans(np.array(v), np.array(t))
#print(v)
#print(vt)
assert(type(vt) == np.ndarray)   # Outer type is preserved
assert(type(vt[0]) == np.double) # Input is ints, but output is np.double
assert(vt.shape == (3, ))


vt = trans(np.array([v]), np.array([t]))
#print(np.array([v]))
#print(vt)
assert(type(vt) == np.ndarray)
assert(vt.shape == (1, 3))
assert(type(vt[0,0]) == np.double)


vt = trans(np.array([v,v]),np.array([t,t]))
#print(np.array([v,v]))
#print(vt)
assert(type(vt) == np.ndarray)
assert(type(vt[0]) == np.ndarray)
assert(type(vt[1]) == np.ndarray)
assert(vt.shape == (2, 3))
for i in range(3):
  assert(type(vt[0,i]) == np.double)
  assert(type(vt[1,i]) == np.double)
#############################################################################


t0 =  [2009,2,2,11,1,1]
t1 =  [[2009,2,2,11,1,1]]
t2 =  [[2009,2,2,11,1,1],[2009,2,2, 1,1,1]]

#############################################################################
v1 = np.array([1,2,3])
vp = trans(v1, t0)
#print(v1)
#print(vp)
assert(vp.shape == v1.shape)

vp = trans(v1, t1)
#print(v1)
#print(vp)
assert(vp.shape == v1.shape)

vp = trans(v1, t2)
#print(v1)
#print(vp)
assert(vp.shape == (2, 3))
#############################################################################

#############################################################################
v1 = np.array([[1,2,3]])

vp = trans(v1, t0)
#print(v1)
#print(vp)
assert(vp.shape == v1.shape)

vp = trans(v1, t1)
#print(v1)
#print(vp)
assert(vp.shape == v1.shape)

vp = trans(v1, t2)
#print(v1)
#print(vp)
assert(vp.shape == (2, 3))
#############################################################################

#############################################################################
v2 = np.array([[1,2,3],[4,5,6]])
vp = trans(v2, t0)
#print(v2)
#print(vp)
assert(vp.shape == (2, 3))

vp = trans(v2, t1)
#print(v2)
#print(vp)
assert(vp.shape == (2, 3))

vp = trans(v2, t2)
#print(v2)
#print(vp)
assert(vp.shape == (2, 3))
#############################################################################
