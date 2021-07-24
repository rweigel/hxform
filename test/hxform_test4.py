import numpy as np
from hxform import hxform as hx

def trans(v,t):
    return hx.transform(v, t, 'GSM','MAG', lib='geopack_08_dp')

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

assert(np.all( trans(v2n  ,t2n)==np.array([trans(vn_0,tn_0),trans(vn_1,tn_1)]) ))

print(trans(v2n  ,t2n))
print(trans(vn_0,tn_0))
print(trans(vn_1,tn_1))
print(trans(v1n_0,tn_0))
print(trans(v1n_1,tn_1))
print(trans(vn_0,t1n_0))
print(trans(vn_1,t1n_1))
print(trans(v1n_0,t1n_0))
print(trans(v1n_1,t1n_1))
