def mlt(pos, time, csys='sph', lib='geopack_08_dp'):
  """Compute magnetic local time given a UT and MAG position or longitude.

  Uses equation 93 in Laundal and Richmond, 2016 (10.1007/s11214-016-0275-y)

  Usage:
  ------
  >>> from hxform import hxform as hx
  >>> mlt = hx.MAGtoMLT(MAGlong, time)
  >>> mlt = hx.MAGtoMLT([MAGlong1, Mlong2, ...], time)

  >>> mlt = hx.MAGtoMLT([MAGx, MAGy, MAGz], time, csys='car')
  >>> mlt = hx.MAGtoMLT([[MAGx1, MAGy1, MAGz1],...], time, csys='car')

  Returns:
  --------
  mlt: float or array-like

  Examples:
  --------
  >>> from hxform import hxform as hx
  >>> mlt = hx.MAGtoMLT(0., [2000, 1, 1, 0, 0, 0])
  >>> print(mlt) # 18.869936573301775

  >>> from hxform import hxform as hx
  >>> mlt = hx.MAGtoMLT([0., 0.], [2000, 1, 1, 0, 0, 0])
  >>> print(mlt) # [18.86993657 18.86993657]

  >>> from hxform import hxform as hx
  >>> mlt = hx.MAGtoMLT([-1., 0., 0.], [2000, 1, 1, 0, 0, 0], csys='car')
  >>> print(mlt) # 6.869936573301775

  >>> from hxform import hxform as hx
  >>> mlt = hx.MAGtoMLT([[-1., 0., 0.],[-1., 0., 0.]], [2000, 1, 1, 0, 0, 0], csys='car')
  >>> print(mlt) # [6.86993657 6.86993657]
  """

  import numpy as np

  from hxform.hxform import transform

  assert csys == 'car' or csys == 'sph', 'csys must be one of ["car", "sph"]'

  pos = np.array(pos)
  time = np.array(time)

  if not isinstance(pos, float):
      pos = np.array(pos)

  if csys == 'sph':
      phi = pos*np.pi/180.
  else:
      if pos.shape == (3, ):
          phi = np.arctan2(pos[1], pos[0])
      else:
          phi = np.arctan2(pos[:, 1], pos[:, 0])

  subsol_pt = transform(np.array([1., 0., 0.]), time, 'GSM', 'MAG', lib=lib)

  if len(subsol_pt.shape) == 1:
      phi_cds = np.arctan2(subsol_pt[1], subsol_pt[0])
  else:
      phi_cds = np.arctan2(subsol_pt[:, 1], subsol_pt[:, 0])

  delta = phi - phi_cds # note np.array([a1, a2, ...])+b == np.array([a1+b, a2+b, ...])

  if isinstance(delta, float):
      delta = np.array([delta])

  idx = np.where(delta > np.pi)
  delta[idx] = delta[idx] - 2.*np.pi
  idx = np.where(delta <= -np.pi)
  delta[idx] = delta[idx] + 2.*np.pi

  if delta.size == 1:
      delta = delta[0]

  MLT = 12. + delta*24./(2.*np.pi)

  return MLT