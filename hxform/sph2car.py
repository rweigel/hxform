def sph2car(*args):
  """Convert r, latitude [degrees], longitude [degrees] to cartesian."""

  import numpy as np
  from utilrsw.np import components2matrix, matrix2components

  try:
    matrix = components2matrix(*args)
    r = matrix[:, 0]
    lat = matrix[:, 1]
    lon = matrix[:, 2]
  except ValueError as e:
    raise e

  if not np.all(r > 0):
    raise ValueError('All radius values must be greater than zero.')

  x = r*np.cos((np.pi/180.0)*lon)*np.cos((np.pi/180.0)*lat)
  y = r*np.sin((np.pi/180.0)*lon)*np.cos((np.pi/180.0)*lat)
  z = r*np.sin((np.pi/180.0)*lat)

  return matrix2components(*args, np.column_stack([x, y, z]))
