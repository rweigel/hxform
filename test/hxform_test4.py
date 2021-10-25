import numpy as np
from hxform import hxform as hx

def trans(v,t):
    return hx.transform(v, t, 'GSM','MAG', lib='geopack_08_dp')

v = [1,2,3]
t = [2009,2,2,11,1,1]

##################################

vt = trans(v,t)
assert(type(vt) == list)
assert(type(vt[0]) == np.double)
assert(len(vt) == 3)

vt = trans(np.array(v),np.array(t))
assert(type(vt) == np.ndarray)
assert(vt.shape == (3,))
assert(type(vt[0]) == np.double)

##################################

##################################
vt = trans([v],[t])
assert(type(vt) == list)
assert(type(vt[0]) == list)
assert(len(vt) == 1)
assert(len(vt[0]) == 3)

vt = trans(np.array([v]),np.array([t]))
assert(type(vt) == np.ndarray)
assert(vt.shape == (1,3))
assert(type(vt[0,0]) == np.double)

##################################

##################################

vt = trans(np.array([v]),np.array([t]))
assert(type(vt) == np.ndarray)
assert(vt.shape == (1,3))
assert(type(vt[0,1]) == np.double)

vt = trans([np.array(v)],[np.array(t)])
assert(type(vt) == list)
assert(len(vt) == 1)
assert(type(vt[0]) == np.ndarray)
assert(type(vt[0][0]) == np.double)

##################################

##################################
vt = trans([v,v],[t,t])
assert(type(vt) == list)
assert(type(vt[0]) == list)
assert(type(vt[1]) == list)
assert(len(vt) == 2)
assert(len(vt[0]) == 3)
assert(len(vt[1]) == 3)
for i in range(3):
    assert(type(vt[0][i]) == np.double)
    assert(type(vt[1][i]) == np.double)

vt = trans(np.array([v,v]),np.array([t,t]))
assert(type(vt) == np.ndarray)
assert(type(vt[0]) == np.ndarray)
assert(type(vt[1]) == np.ndarray)
assert(vt.shape == (2, 3))
for i in range(3):
    assert(type(vt[0,i]) == np.double)
    assert(type(vt[1,i]) == np.double)

##################################
