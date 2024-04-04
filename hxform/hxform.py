def transform(v, time, csys_in, csys_out, ctype_in='car', ctype_out='car', lib='cxform'):
  """Transform between coordinates systems using cxform, Geopack, SpacePy, SpiceyPy, SSCWeb, or SunPy.

  Parameters
  ----------
  v : array-like

      (Nv, 3) float np.array

      np.array of three floats

      list of three floats

      list containing lists of three floats

      list of 3-element np.arrays

  time : array-like
          list of 3+ ints
          list containing lists of 3+ ints
          np.array of 3+ ints
          (Nt, 3) float np.array, where Nt = 1 or Nt = Nv

          The 3+ ints are [year, month, day, [hours, [minutes, [seconds]]]]
          Zeros are used for any missing optional value.

  csys_in : str
            One of MAG, GEI, GEO, GSE, GSM, SM

  csys_out : str
              One of MAG, GEI, GEO, GSE, GSM, SM

  ctype_in : str
              'car' (default) or 'sph'
              For spherical coordinates, `v` should be in r, latitude, longitude,
              with angles in degrees.

  ctype_out : str
              'car' (default) or 'sph'

  lib : str
        'cxform' (default), 'geopack_08_dp', 'spacepy', or 'spacepy-irbem'

  Returns
  -------
  array-like with dimensions matching either `time` (if `Nt` != 1 and `Nv` = 1) or
  `v` (if `Nv` =! 1 and `Nt` = 1). If `Nv` and `Nt` != 1, dimensions are same as `v`.

  Return type will match that of `v`. Note that if a list of 3-element np.arrays are
  passed, execution time will be long. Use `np.ndarrays` for `v` and `time` for fastest
  execution time.

  Examples
  --------
  >>> from hxform import hxform as hx
  >>> t1 = [2000, 1, 1, 0, 0, 0] # or np.array([2000, 1, 1, 0, 0, 0])
  >>> v1 = [0, 0, 1]             # or np.array([0, 0, 1])
  >>> # All of the following are equivalent and return a list with three floats
  >>> from hxform import hxform as hx
  >>> hx.transform(v1, time1, 'GSM', 'GSE')
  >>> hx.transform(v1, time1, 'GSM', 'GSE', ctype_in='car')
  >>> hx.transform(v1, time1, 'GSM', 'GSE', ctype_in='car', ctype_out='car')

  The following 3 calls return a list with two lists of 3 elements

  1. Transform two vectors at same time t1

      >>> from hxform import hxform as hx
      >>> hx.transform([v1, v1], t1, 'GSM', 'GSE')

  2. Transform one vector at two different times

      >>> from hxform import hxform as hx
      >>> hx.transform(v1, [t1, t1], 'GSM', 'GSE')

  3. Transform two vectors, each at different times

      >>> from hxform import hxform as hx
      >>> hx.transform([v1, v1], [t1, t1], 'GSM', 'GSE')
  """
  import numpy as np
  import hxform

  assert lib in hxform.info.known_libs(), 'lib must be one of {}'.format(hxform.info.known_libs())
  assert ctype_in in ['car', 'sph'], 'ctype_in must be one of ["car", "sph"]'
  assert ctype_out in ['car', 'sph'], 'ctype_out must be one of ["car", "sph"]'

  info = hxform.info.lib_info(lib)
  assert csys_in in info['systems'], 'For lib={}, csys_in must be one of {}'.format(lib, info['systems'])
  assert csys_out in info['systems'], 'For lib={}, csys_out must be one of {}'.format(lib, info['systems'])

  v_outertype = type(v)
  v_innertype = type(v[0])

  v = np.array(v, dtype=np.double)
  try:
    time = np.array(time, dtype=np.int32)
  except:
    # This will occur for input of the form [[2000, 1, 1], [2000, 1, 1, 1]]
    # where not all time values have the same number of elements.
    for i in range(len(time)):
      print(hxform.timelib.tpad(time[i], length=6))
      time[i] = hxform.timelib.tpad(time[i], length=6)
    time = np.array(time, dtype=np.int32)

  if len(time.shape) == 1:
    time = np.array([time])
  if len(v.shape) == 1:
    v = np.array([v])

  Nv = v.shape[0]     # Number of vectors
  Nt = time.shape[0]  # Number of times

  assert(len(v[0]) == 3)
  assert(len(time.shape) == 2 and len(v.shape) == 2)
  assert(Nv == Nt or Nt == 1 or Nv == 1)

  if csys_in == csys_out or lib.startswith('spacepy') or lib.startswith('spiceypy') or lib == 'sscweb' or lib == 'sunpy':
    # TODO: Could avoid expanding time or v by putting logic to test for
    # Nt == 1 or Nv == 1 in transform loops for each lib.
    from numpy import matlib
    if Nt == 1 and Nv > 1:
      time = matlib.repmat(time, Nv, 1)

    if Nv == 1 and Nt > 1:
      v = matlib.repmat(v, Nt, 1)

  # vp means "v prime", which is vector v in cys_in transformed to csys_out
  if lib.startswith('spiceypy') or lib == 'sscweb' or lib == 'sscweb' or lib == 'sunpy':
    vp = np.full(v.shape, np.nan)

  if ctype_in == 'sph' and lib != 'sscweb':
    v[:,0], v[:,1], v[:,2] = StoC(v[:,0], v[:,1], v[:,2])

  import time as time_
  execution_start = time_.time()

  if lib == 'cxform':
    import os
    import glob
    import ctypes
    this_script = os.path.join(os.path.dirname(__file__), "..", "hxform", "cxform_wrapper*")

    for lib_file in glob.glob(this_script):
      # The name of the .so or .dll file will not be the same on all
      # systems, so we need to find it. (For example, on one system
      # it is cxform_wrapper.cpython-37m-darwin.so.)
      # TODO: Find a better way to do this.
      break
    lib_path = os.path.join(lib_file)
    lib_obj = ctypes.cdll.LoadLibrary(lib_path)

    Nt = time.shape[0]

    if time.shape[1] < 3:
      raise ValueError("At least year, month, and day must be given for time.")

    nz = time.shape[1]
    if nz != 6:
      # Pad time. TODO: Do this in wrapper so extra memory is not needed.
      tmp = np.zeros((time.shape[0], 6-nz), dtype=np.int32)
      time = np.concatenate((time, tmp), 1)

    if Nt == 1:
      vp = np.full(v.shape, np.nan)
    else:
      vp = np.full((Nt, 3), np.nan)

    ret = lib_obj.cxform_wrapper(
            ctypes.c_void_p(v.ctypes.data),
            ctypes.c_void_p(time.ctypes.data),
            ctypes.c_char_p(str.encode(csys_in)),
            ctypes.c_char_p(str.encode(csys_out)),
            ctypes.c_void_p(vp.ctypes.data),
            ctypes.c_int(v.shape[0]),
            ctypes.c_int(int(Nt))
        )

  if lib == 'geopack_08_dp':

    import hxform.geopack_08_dp as geopack_08_dp
    trans = csys_in + 'to' + csys_out
    dtime = np.array(hxform.timelib.ints2doy(time))

    if v.shape[0] <= time.shape[0]:
      outsize = time.shape[0]
    else:
      outsize = v.shape[0]

    if ctype_in == 'sph':
      v[:,0], v[:,1], v[:,2] = StoC(v[:,0], v[:,1], v[:,2])

    vp = geopack_08_dp.transform(v, trans, dtime, outsize)

    if ctype_out == 'sph':
      vp[:,0], vp[:,1], vp[:,2] = CtoS(vp[:,0], vp[:,1], vp[:,2])

  if lib.startswith('spacepy'):
    import spacepy.coordinates as sc
    from spacepy.time import Ticktock

    if lib.endswith('-irbem'):
      cvals = sc.Coords(v, csys_in, ctype_in, use_irbem=True)
    else:
      cvals = sc.Coords(v, csys_in, ctype_in, use_irbem=False)

    if len(time.shape) == 1:
      # SpacePy requires time values to be strings with 1-second precision
      t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.tpad(time, length=6))
    else:
      t_str = []
      for i in range(time.shape[0]):
        t_str.append('%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.tpad(time[i,:], length=6)))
      t_str = np.array(t_str)

    cvals.ticks = Ticktock(t_str, 'ISO')
    vp = cvals.convert(csys_out, ctype_out).data

  if lib == 'sunpy':

    import astropy.coordinates
    import sunpy.coordinates

    if ctype_in == 'sph':
      units = [astropy.units.m, astropy.units.deg, astropy.units.deg]
      representation_type = 'spherical'
    else:
      representation_type = 'cartesian'
      one = astropy.constants.R_earth
      # TODO: Use
      #   one = astropy.units.m
      # when https://github.com/sunpy/sunpy/pull/7530 is merged.
      units = [one, one, one]

    frame_in = info['system_aliases'][csys_in]
    frame_out = info['system_aliases'][csys_out]

    if Nt == 1:
      obstime = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.tpad(time[0,:], length=6))
      kwargs = {
        "x": v[:,0]*units[0],
        "y": v[:,1]*units[1],
        "z": v[:,2]*units[2],
        "frame": frame_in,
        "obstime": obstime,
        "representation_type": representation_type
      }
      coord = astropy.coordinates.SkyCoord(**kwargs)

      if ctype_out == 'car':
        coord = coord.transform_to(frame_out).cartesian/one
      if ctype_out == 'sph':
        coord = coord.transform_to(frame_out).spherical/one

      vp = coord.xyz.decompose().value.transpose()

    else:
      for i in range(Nt):
        obstime = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.tpad(time[i,:], length=6))
        kwargs = {
          "x": v[i,0]*units[0],
          "y": v[i,1]*units[1],
          "z": v[i,2]*units[2],
          "frame": frame_in,
          "obstime": obstime,
          "representation_type": representation_type
        }

        coord_in = astropy.coordinates.SkyCoord(**kwargs)

        if ctype_out == 'car':
          coord_out = coord_in.transform_to(frame_out).cartesian/one
        if ctype_out == 'sph':
          coord_out = coord_in.transform_to(frame_out).spherical/one

        vp[i,:] = coord_out.xyz.decompose().value

  if lib.startswith('spiceypy'):

    import os
    import spiceypy

    kernels = info['kernels']
    for kernel in kernels:
      rel_path = os.path.join('..','demo','spiceypy', 'kernels', kernel)
      kernel_file = os.path.join(os.path.dirname(__file__), rel_path)
      spiceypy.furnsh(kernel_file)

    for i in range(time.shape[0]):
      time_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.tpad(time[i,:], length=6))
      et = spiceypy.str2et(time_str)
      matrix = spiceypy.pxform(csys_in, csys_out, et)
      vp[i,:] = spiceypy.mxv(matrix, v[i,:])

    spiceypy.kclear()

  if lib == 'sscweb':
    import re
    import requests
    import time as time_

    if csys_in in ['GEO', 'GM'] and ctype_in == 'car':
      v[:,0], v[:,1], v[:,2] = CtoS(v[:,0], v[:,1], v[:,2])

    for i in range(time.shape[0]):

      time_.sleep(0.1)
      year = time[i][0]
      doy_ = hxform.timelib.doy(time[i][0:3])
      time_str = f'{year} {doy_:03d} {time[i][3]:02d}:{time[i][4]:02d}:{time[i][5]:02d}'
      url = "https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi?"
      if csys_in in ['GEO', 'GM']:
        url += f"epoch={time_str}&x=&y=&z=&lat={v[i][1]:f}&lon={v[i][2]:f}&r={v[i][0]:f}&action={csys_in}"
      else:
        url += f"epoch={time_str}&x={v[i][0]:f}&y={v[i][1]:f}&z={v[i][2]:f}&lat=&lon=&r=&action={csys_in}"

      # This is likely to fail if many requests are made.
      try:
        r = requests.get(url)
      except:
        raise Exception(f"Failed to fetch URL: {url}.")
        return None

      if r.status_code != 200:
        raise Exception(f"Failed to fetch URL: {url}. Status code: {r.status_code}.")

      text = r.text

      # Split the text into lines
      lines = text.split("\n")

      start = None
      end = None
      for idx, line in enumerate(lines):
        if re.search('^ Radial distance', line):
          start = idx
        if re.search('^ REGION', line):
          end = idx - 1
          break

      if start is None:
        raise Exception(f"Failed to find start of table in URL: {url}. Returned HTML:\n{text}")

      # Extract the table lines
      R = float(lines[start].split()[2].strip())
      table_lines = lines[start+2:end]

      # Parse the table into a list of dictionaries
      result = {}
      for line in table_lines:
        parts = line.split()
        result[parts[0]] = {
            "R": R,
            "Lat": float(parts[1]),
            "Long": float(parts[2]),
            "X": float(parts[3]),
            "Y": float(parts[4]),
            "Z": float(parts[5]),
        }
        if len(parts) > 6:
          result[parts[0]]["hh.hhhhh"] = float(parts[6])

      if ctype_out == 'car':
        vp[i][0] = result[csys_out]['X']
        vp[i][1] = result[csys_out]['Y']
        vp[i][2] = result[csys_out]['Z']
      else:
        vp[i][0] = result[csys_out]['R']
        vp[i][1] = result[csys_out]['Lat']
        vp[i][2] = result[csys_out]['Long']

  if ctype_out == 'sph' and lib != 'sscweb':
    vp[:,0], vp[:,1], vp[:,2] = CtoS(vp[:,0], vp[:,1], vp[:,2])

  transform.execution_time = time_.time() - execution_start

  if issubclass(v_outertype, np.ndarray):
    if Nv == 1 and Nt == 1 and not issubclass(v_innertype, np.ndarray):
      return vp[0]
    return vp
  elif issubclass(v_innertype, np.ndarray):
    return v_outertype(vp)
  else:
    if Nv == 1 and Nt == 1 and v_innertype != list:
      return vp.tolist()[0]
    return vp.tolist()


def transform_matrix(time, csys_in, csys_out, lib='geopack_08_dp'):
  import numpy as np
  kwargs = {'ctype_in': 'car', ctype_out: 'car', 'lib': lib}
  b1 = transform(np.array([1., 0., 0.]), time, csys_in, csys_out, **kwargs)
  b2 = transform(np.array([0., 1., 0.]), time, csys_in, csys_out, **kwargs)
  b3 = transform(np.array([0., 0., 1.]), time, csys_in, csys_out, **kwargs)
  return np.column_stack([b1, b2, b3])


def CtoS(x, y, z):
    """Convert from cartesian to r, latitude [degrees], longitude [degrees]."""
    import numpy as np
    r = np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
    assert np.all(r > 0), 'radius must be greater than zero.'
    lat = 90.0 - (180.0/np.pi)*np.arccos(z/r)
    lon = (180.0/np.pi)*np.arctan2(y, x)

    return r, lat, lon


def StoC(r, lat, lon):
    """Convert r, latitude [degrees], longitude [degrees] to cartesian."""
    import numpy as np
    assert np.all(r > 0), 'radius must be greater than zero.'
    x = r*np.cos((np.pi/180.0)*lon)*np.cos((np.pi/180.0)*lat)
    y = r*np.sin((np.pi/180.0)*lon)*np.cos((np.pi/180.0)*lat)
    z = r*np.sin((np.pi/180.0)*lat)

    return x, y, z


def get_spherical_vector_components(v_cart, x_cart):

  import numpy as np
  x = x_cart[:,0]
  y = x_cart[:,1]
  z = x_cart[:,2]
  vx = v_cart[:,0]
  vy = v_cart[:,1]
  vz = v_cart[:,2]

  r = np.sqrt(x**2 + y**2 + z**2)
  L = np.sqrt(x**2 + y**2)

  v_r = (x*vx + y*vy + z*vz)/r
  v_theta = ((x*vx + y*vy)*z - (L**2)*vz)/(r*L)
  v_phi = (-y*vx + x*vy)/L

  return np.column_stack([v_r, v_theta, v_phi])


def get_NED_vector_components(v_cart, x_cart):
  import numpy as np
  v_sph = get_spherical_vector_components(v_cart, x_cart)
  v_north = -v_sph[:,1]
  v_east  = +v_sph[:,2]
  v_down  = -v_sph[:,0]
  return np.column_stack([v_north, v_east, v_down])


def MAGtoMLT(pos, time, csys='sph', lib='geopack_08_dp'):
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
