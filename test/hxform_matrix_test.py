import numpy
import hxform

t1 = numpy.array([2010, 12, 30, 0, 0, 0])
v1 = numpy.array([1., 1., 1.])
v2 = numpy.array([v1, 2*v1])
t2 = numpy.array([t1, t1])
frame_in = 'MAG'
frame_out = 'GSE'
kwargs = {
    'frame_in': frame_in,
    'frame_out': frame_out,
}

#vt = hxform.transform(v1, t1, **kwargs)
#matrix = hxform.matrix(t1, **kwargs)
#assert (numpy.all(numpy.abs(vt == numpy.dot(matrix, v1)) < 1e-15))

libs = hxform.libs()
for lib in libs:
  kwargs['lib'] = lib
  if lib == 'sunpy':
    continue
  print(lib)
  # Transform single vector at single timestamp
  vt = hxform.transform(v1, t1, **kwargs)
  matrix = hxform.matrix(t1, **kwargs)
  print(matrix)
  diff = vt - numpy.dot(matrix, v1)
  print(diff)
  assert numpy.all(numpy.abs(diff) < 1e-15)

  continue
  # Transform two vectors at with same timestamp
  vt = hxform.transform(v2, t1, **kwargs)
  matrix = hxform.matrix(t1, **kwargs)
  assert numpy.all(numpy.abs(vt[0, :] - numpy.dot(matrix, v2[0, :])) < 1e-15)
  assert numpy.all(numpy.abs(vt[1, :] - numpy.dot(matrix, v2[1, :])) < 1e-15)

  # Transform two vectors at two timestamps
  vt = hxform.transform(v2, t2, **kwargs)
  matrices = hxform.matrix(t2, **kwargs)
  for i in range(matrices.shape[0]):
    assert (numpy.all(numpy.abs(vt[i, :] - numpy.dot(matrices[i, :, :], v2[i, :])) < 1e-15))