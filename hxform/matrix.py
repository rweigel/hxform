def matrix(t, frame_in, frame_out, lib='geopack_08_dp'):
  import numpy

  from hxform.hxform import transform

  kwargs = {'ctype_in': 'car', 'ctype_out': 'car', 'lib': lib}
  c1 = transform(numpy.array([1., 0., 0.]), t, frame_in, frame_out, **kwargs)
  c2 = transform(numpy.array([0., 1., 0.]), t, frame_in, frame_out, **kwargs)
  c3 = transform(numpy.array([0., 0., 1.]), t, frame_in, frame_out, **kwargs)
  if len(c1.shape) == 1:
    m = numpy.column_stack([c1, c2, c3])
    if isinstance(t, list):
      return m.tolist()
    return m

  m = numpy.full((t.shape[0], 3, 3), numpy.nan)
  for i in range(t.shape[0]):
    m[i, :, 0] = c1[i, :]
    m[i, :, 1] = c2[i, :]
    m[i, :, 2] = c3[i, :]
  if isinstance(t, list):
    return m.tolist()

  return m
