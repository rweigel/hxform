# API Test

import numpy as np
from hxform import hxform as hx

def trans(v,t):
    return hx.transform(v, t, 'GSM','MAG', lib='geopack_08_dp')

v = [1,2,3]
t = [2009,2,2,11,1,1]
out = hx.transform(v, t, 'GSM','MAG', lib='geopack_08_dp')
##################################

vt = trans(v,t)
assert(type(vt) == list)
assert(type(vt[0]) == np.double)
assert(len(vt) == 3)

vt = trans(np.array(v),np.array(t))
assert(type(vt) == np.ndarray)
assert(type(vt[0]) == np.double)
assert(vt.shape == (3,))

##################################

##################################
vt = trans([v],[t])
#print(vt)
#print(type(vt))
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

vt = trans(np.array(v),np.array(t))
assert(type(vt) == np.ndarray)
assert(vt.shape == (3,))
for i in range(3):
    assert(type(vt[i]) == np.double)


vt = trans(np.array([v]),np.array([t]))
assert(type(vt) == np.ndarray)
assert(vt.shape == (1,3))

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

v2n   = np.array([[1,2,3],[4,5,6]])
v1n_0 = np.array([[1,2,3]])
v1n_1 = np.array([[4,5,6]])
vn_0  = np.array([1,2,3])
vn_1  = np.array([4,5,6])

t2n   =  [[2009,2,2,11,1,1],[2009,2,2, 1,1,1]]
t1n_0 =  [[2009,2,2,11,1,1]]
t1n_1 =  [[2009,2,2, 1,1,1]]
tn_0  =  [2009,2,2,11,1,1]
tn_1  =  [2009,2,2, 1,1,1]

assert( trans(v2n  ,t2n).shape == (2,3) )
assert( trans(v1n_0,t2n).shape == (2,3) )
assert( trans(vn_0 ,t2n).shape == (2,3) )

assert( trans(v2n  ,t1n_0).shape == (2,3) )
assert( trans(v1n_0,t1n_0).shape == (1,3) )
assert( trans(vn_0 ,t1n_0).shape == (1,3) )

assert( trans(v2n  ,tn_0).shape == (2,3) )
assert( trans(v1n_0,tn_0).shape == (1,3) )
assert( trans(vn_0 ,tn_0).shape == ( 3,) )

assert(np.all( trans(v2n  ,t2n) == np.array([trans(vn_0,tn_0),trans(vn_1,tn_1)]) ))
assert(np.all( trans(vn_0,tn_0) == trans(v1n_0,t1n_0).ravel() ))
assert(np.all( trans(vn_0,tn_0) == trans(v1n_0, tn_0).ravel() ))
assert(np.all( trans(vn_0,tn_0) == trans( vn_0,t1n_0).ravel() ))
assert(np.all( trans(vn_1,tn_1) == trans(v1n_1,t1n_1).ravel() ))
assert(np.all( trans(vn_1,tn_1) == trans(v1n_1, tn_1).ravel() ))
assert(np.all( trans(vn_1,tn_1) == trans( vn_1,t1n_1).ravel() ))
