import numpy
import hxform

t1  = numpy.array([2010, 12, 30, 0, 0, 0])
v1 = numpy.array([1., 1., 1.])
v2 = numpy.array([v1, 2*v1])
t2 = numpy.array([t1, t1])
frame_in = 'GEO'
frame_out = 'GSE'
kwargs = {
    'frame_in': frame_in,
    'frame_out': frame_out,
    'lib': 'geopack_08_dp'
}

vt = hxform.transform(v1, t1, **kwargs)
matrix = hxform.matrix(t1, **kwargs)
assert (numpy.all(numpy.abs(vt == numpy.dot(matrix, v1)) < 1e-15))

# Transform two vectors at with same timestamp
vt = hxform.transform(v2, t1, **kwargs)
matrix = hxform.matrix(t1, **kwargs)
assert (numpy.all(numpy.abs(vt == numpy.dot(matrix, v2[0, :])) < 1e-15))
assert (numpy.all(numpy.abs(vt == numpy.dot(matrix, v2[1, :])) < 1e-15))

# Transform two vectors at two timestamps
vt = hxform.transform(v2, t2, **kwargs)
matrices = hxform.matrix(t2, **kwargs)
for i in range(matrices.shape[0]):
  assert (numpy.all(numpy.abs(vt[i, :] == numpy.dot(matrices[i, :, :], v2[i, :])) < 1e-15))