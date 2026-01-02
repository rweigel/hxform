def mlt(pos, time, csys='sph', lib='geopack_08_dp'):
  """Compute magnetic local time given a UT and MAG position or longitude.

  MLT [fractional hours] = 12 + (phi - phi_s)*24/360

  where phi_s is the longitude of the subsolar point in the centered dipole
  frame and phi is position in the centered dipole frame. Both angles are given
  in degrees. The MAG frame is used for the centered dipole frame.

  Usage:
  ------
  hxform.mlt(MAGlong, time)
  hxform.mlt([MAGlong1, Mlong2, ...], time)
  hxform.mlt([MAGx, MAGy, MAGz], time, csys='car')
  hxform.mlt([[MAGx1, MAGy1, MAGz1],...], time, csys='car')

  Returns:
  --------
  mlt: float or array-like

  Examples:
  --------
  >>> import hxform
  >>> mlt = hxform.mlt(0., [2000, 1, 1, 0, 0, 0])
  >>> print(mlt) # 18.8708344

  >>> mlt = hxform.mlt([0., 0.], [2000, 1, 1, 0, 0, 0])
  >>> print(mlt) # [18.8708344 18.8708344]

  >>> mlt = hxform.mlt([-1., 0., 0.], [2000, 1, 1, 0, 0, 0], csys='car')
  >>> print(mlt) # 6.869936573301775

  >>> mlt = hxform.mlt([[-1., 0., 0.],[-1., 0., 0.]], [2000, 1, 1, 0, 0, 0], csys='car')
  >>> print(mlt) # [6.86993657 6.86993657]
  """

  import hxform

  import numpy as np

  if csys not in ['car', 'sph']:
    raise ValueError(f'csys must be one of ["car", "sph"], got {csys}')

  if lib not in hxform.libs():
    raise ValueError(f"Library '{lib}' not in available libraries: {hxform.libs()}")

  pos = np.array(pos)

  if not isinstance(pos, float):
    pos = np.array(pos)

  if csys == 'sph':
    phi = pos*np.pi/180.
  else:
    if pos.shape == (3, ):
      phi = np.arctan2(pos[1], pos[0])
    else:
      phi = np.arctan2(pos[:, 1], pos[:, 0])

  subsol_pt = hxform.subsolar_point(np.array(time), lib=lib)

  if len(subsol_pt.shape) == 1:
    phi_cds = np.arctan2(subsol_pt[1], subsol_pt[0])
  else:
    phi_cds = np.arctan2(subsol_pt[:, 1], subsol_pt[:, 0])

  delta = phi - phi_cds

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


if __name__ == '__main__':

  import hxform
  lib = 'geopack_08_dp'
  t = '2000-01-01T00:00:00'

  mlong = 0.0
  mlt = hxform.mlt(mlong, t, lib=lib)
  print(mlt) # 18.8708344

  mlt = hxform.mlt([mlong, mlong], t, lib=lib)
  print(mlt) # [18.8708344 18.8708344]

  mlt = hxform.mlt([-1., 0., 0.], t, lib=lib, csys='car')
  print(mlt) # 6.8708344027652934

  mlt = hxform.mlt([[-1., 0., 0.], [-1., 0., 0.]], [2000, 1, 1, 0, 0, 0], lib=lib, csys='car')
  print(mlt) # [6.8708344 6.8708344]

  libs = hxform.libs()
  print(f'MLT at 0 deg MAG long at {t}Z')
  for lib in libs:
    try:
      mlt = hxform.mlt(0., [2000, 1, 1, 0, 0, 0], lib=lib)
      print(f'{lib:13s}: {mlt}')
    except Exception as e:
      print(f'{lib:13s}: Error: {e}')
