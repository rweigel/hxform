def car2sph(*args):
  """Convert from cartesian to r, latitude [degrees], longitude [degrees]."""

  import numpy as np
  from utilrsw.np import components2matrix, matrix2components

  try:
    matrix = components2matrix(*args)
    x = matrix[:, 0]
    y = matrix[:, 1]
    z = matrix[:, 2]
  except ValueError as e:
    raise e

  r = np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
  assert np.all(r > 0), 'radius must be greater than zero.'
  lat = 90.0 - (180.0/np.pi)*np.arccos(z/r)
  lon = (180.0/np.pi)*np.arctan2(y, x)

  return matrix2components(*args, np.column_stack([r, lat, lon]))
