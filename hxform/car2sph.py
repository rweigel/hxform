def car2sph(*args):
  """Convert from cartesian to r, latitude [degrees], longitude [degrees].
  Usage:
    r, lat, lon = car2sph(x, y, z) where x, y, z are 1D arrays, lists, or
    tuples of the same length.

    rtp = car2sph(np.array(xyz)) where xyz is a 2D array-like object with
    shape (N, 3). xyz may be a list/tuple of lists/tuples or 1D numpy arrays.
    rtp (r, lat, lon) has the same structure as the input.
  """

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
  if not np.all(r > 0):
    raise ValueError('x^2 + y^2 + z^2 must be greater than zero.')

  lat = 90.0 - (180.0/np.pi)*np.arccos(z/r)
  lon = (180.0/np.pi)*np.arctan2(y, x)

  return matrix2components(*args, np.column_stack([r, lat, lon]))
