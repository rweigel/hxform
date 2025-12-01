def transform(v, t, frame_in, frame_out, ctype_in='car', ctype_out='car', lib='cxform'):
  """Transform between coordinates frames using cxform, Geopack-08, PySPEDAS, 
  SpacePy, SpiceyPy, SSCWeb, or SunPy.

  Parameters
  ----------
  v : array-like

      (Nv, 3) float np.array

      np.array of three floats

      list of three floats

      list containing lists of three floats

      list of 3-element np.arrays

  t : array-like

      list of 3+ ints ([year, month, day, [hours, [minutes, [seconds]]]])
      Zeros are used for any missing optional value.

      (Nt, 3+) float np.array, where Nt = 1 or Nt = Nv

      list containing lists of 3+ ints

      np.array of 3+ ints

  lib : str
        Library to use for the transformation. See hxform.libs() for list.

  frame_in : str
             See hxform.frames(lib)

  frame_out : str
              One of MAG, GEI, GEO, GSE, GSM, SM

  ctype_in : str
            'car' (default) or 'sph'
            For spherical coordinates, `v` is r, latitude, longitude, with
            angles in degrees.

  ctype_out : str
              'car' (default) or 'sph'


  Returns
  -------
  array-like with dimensions matching either `t` (if `Nt` != 1 and `Nv` = 1) or
  `v` (if `Nv` =! 1 and `Nt` = 1). If `Nv` and `Nt` != 1, dimensions are same as `v`.

  Return type will match that of `v`. Note that if a list of 3-element np.arrays are
  passed, execution time will be long. Use `np.ndarrays` for `v` and `t` for fastest
  execution time.

  Examples
  --------
  >>> from hxform import hxform as hx
  >>> t1 = [2000, 1, 1, 0, 0, 0] # or np.array([2000, 1, 1, 0, 0, 0])
  >>> v1 = [0, 0, 1]             # or np.array([0, 0, 1])
  >>> # All of the following are equivalent and return a list with three floats
  >>> from hxform import hxform as hx
  >>> hx.transform(v1, t1, 'GSM', 'GSE')
  >>> hx.transform(v1, t1, 'GSM', 'GSE', ctype_in='car')
  >>> hx.transform(v1, t1, 'GSM', 'GSE', ctype_in='car', ctype_out='car')

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

  assert lib in hxform.libs(), f'lib must be one of {hxform.info.libs()}'
  assert ctype_in in ['car', 'sph'], 'ctype_in must be one of ["car", "sph"]'
  assert ctype_out in ['car', 'sph'], 'ctype_out must be one of ["car", "sph"]'

  info = hxform.lib_info(lib)
  for csys in [frame_in, frame_out]:
    emsg = f'For lib={lib}, {csys} must be one of {info["frames"]}'
    assert frame_in in info['frames'], emsg

  v_outertype = type(v)
  v_innertype = type(v[0])

  v = np.array(v, dtype=np.double)
  try:
    t = np.array(t, dtype=np.int32)
  except:
    # This will occur for input of the form [[2000, 1, 1], [2000, 1, 1, 1]]
    # where not all time values have the same number of elements.
    for i in range(len(t)):
      t[i] = hxform.timelib.ints_pad(t[i], length=6)
    t = np.array(t, dtype=np.int32)

  if len(t.shape) == 1:
    t = np.array([t])
  if len(v.shape) == 1:
    v = np.array([v])

  if t.shape[1] > 6:
    # Keep only year, month, day, hour, minute, second (drop microseconds)
    t = t[:,0:6]

  Nv = v.shape[0]  # Number of vectors
  Nt = t.shape[0]  # Number of times

  assert(len(v[0]) == 3)
  assert(len(t.shape) == 2 and len(v.shape) == 2)
  assert(Nv == Nt or Nt == 1 or Nv == 1)

  tile = False
  tile = frame_in == frame_out
  tile = tile or lib.startswith('spacepy') or lib.startswith('spiceypy')
  tile = tile or lib in ['sscweb', 'sunpy', 'pyspedas']
  if tile:
    if Nt == 1 and Nv > 1 and lib != 'sunpy':
      # In this case, we could compute matrix once for the single time and
      # apply it to all times. This would speed up execution time for some libs.
      # However, many failures in computing the matrix for certain libraries;
      # see hxform/test/matrix_test.py
      t = np.tile(t, (Nv, 1))

    if Nv == 1 and Nt > 1:
      v = np.tile(v, (Nt, 1))

  if ctype_in == 'sph':
    v[:,0], v[:,1], v[:,2] = sph2car(v[:,0], v[:,1], v[:,2])


  if lib == 'cxform':
    vt, execution_time = _cxform(v, t, frame_in, frame_out)


  if lib == 'geopack_08_dp':
    vt, execution_time = _geopack_08_dp(v, t, frame_in, frame_out)


  if lib == 'pyspedas':
    vt, execution_time = _pyspedas(v, t, frame_in, frame_out)


  if lib.startswith('spacepy'):
    vt, execution_time = _spacepy(v, t, frame_in, frame_out, lib)


  if lib == 'sunpy':
    vt, execution_time = _sunpy(v, t, frame_in, frame_out)


  if lib.startswith('spiceypy'):
    vt, execution_time = _spiceypy(v, t, frame_in, frame_out, lib)


  if lib == 'sscweb':
    vt, execution_time = _sscweb(v, t, frame_in, frame_out)

  if False and Nt == 1 and Nv > 1:
    matrix = hxform.matrix(t[0,:], frame_in=frame_in, frame_out=frame_out, lib=lib)
    vt_alt = np.dot(v, matrix)
    import pdb; pdb.set_trace()
    print(np.max(np.abs(vt - vt_alt)))

  transform.execution_time = execution_time

  if ctype_out == 'sph':
    vt[:,0], vt[:,1], vt[:,2] = car2sph(vt[:,0], vt[:,1], vt[:,2])

  if issubclass(v_outertype, np.ndarray):
    if Nv == 1 and Nt == 1 and not issubclass(v_innertype, np.ndarray):
      return vt[0]
    return vt
  elif issubclass(v_innertype, np.ndarray):
    return v_outertype(vt)
  else:
    if Nv == 1 and Nt == 1 and v_innertype is not list:
      return vt.tolist()[0]
    return vt.tolist()


def _cxform(v, t, frame_in, frame_out):
  import os
  import glob
  import ctypes
  import time
  import numpy as np

  this_script = os.path.join(os.path.dirname(__file__), "..", "hxform", "cxform_wrapper*")

  for lib_file in glob.glob(this_script):
    # The name of the .so or .dll file will not be the same on all
    # frames, so we need to find it. (For example, on one system
    # it is cxform_wrapper.cpython-37m-darwin.so.)
    # TODO: Find a better way to do this.
    break
  lib_path = os.path.join(lib_file)
  lib_obj = ctypes.cdll.LoadLibrary(lib_path)

  Nt = t.shape[0]

  if t.shape[1] < 3:
    raise ValueError("At least year, month, and day must be given for time.")

  execution_start = time.time()
  nz = t.shape[1]
  if nz != 6:
    # Pad time. TODO: Do this in wrapper so extra memory is not needed.
    tmp = np.zeros((t.shape[0], 6-nz), dtype=np.int32)
    t = np.concatenate((t, tmp), 1)

  if Nt == 1:
    vt = np.full(v.shape, np.nan)
  else:
    vt = np.full((Nt, 3), np.nan)

  ret = lib_obj.cxform_wrapper(
          ctypes.c_void_p(v.ctypes.data),
          ctypes.c_void_p(t.ctypes.data),
          ctypes.c_char_p(str.encode(frame_in)),
          ctypes.c_char_p(str.encode(frame_out)),
          ctypes.c_void_p(vt.ctypes.data),
          ctypes.c_int(v.shape[0]),
          ctypes.c_int(int(Nt))
      )
  execution_stop = time.time()

  return vt, execution_stop - execution_start


def _geopack_08_dp(v, t, frame_in, frame_out):
  import time
  import numpy as np
  import hxform
  import hxform.geopack_08_dp as geopack_08_dp

  trans = frame_in + 'to' + frame_out
  dtime = np.array(hxform.timelib.ints2doy(t))

  if v.shape[0] <= t.shape[0]:
    outsize = t.shape[0]
  else:
    outsize = v.shape[0]

  execution_start = time.time()
  vt = geopack_08_dp.transform(v, trans, dtime, outsize)
  execution_stop = time.time()

  return vt, execution_stop - execution_start


def _pyspedas(v, t, frame_in, frame_out):
  import os
  import time
  import numpy as np
  import hxform

  os.environ["PYSPEDAS_LOGGING_LEVEL"] = "error"
  from pyspedas import cotrans
  from pyspedas import time_double

  time_in = []
  for i in range(t.shape[0]):
    tstr = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.ints_pad(t[i,:], length=6))
    time_in.append(time_double(tstr))

  execution_start = time.time()
  vt = cotrans(time_in=time_in, data_in=v, coord_in=frame_in, coord_out=frame_out)
  if isinstance(vt, list):
    # https://github.com/spedas/pyspedas/issues/1273
    vt = np.array(vt)
  execution_stop = time.time()

  return vt, execution_stop - execution_start


def _spiceypy(v, t, frame_in, frame_out, lib):

  import os
  import time

  import numpy as np
  import spiceypy
  import hxform

  vt = np.full(v.shape, np.nan)

  info = hxform.lib_info(lib)
  kernel_dir = info['kernel']['dir']
  kernel_files = info['kernel']['files']
  for kernel in kernel_files:
    kernel_file = os.path.join(kernel_dir, kernel)
    if not os.path.isfile(kernel_file):
      raise FileNotFoundError(f"Required SPICE kernel file not found: {kernel_file}")
    spiceypy.furnsh(kernel_file)

  execution_start = time.time()
  for i in range(t.shape[0]):
    time_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.ints_pad(t[i,:], length=6))
    et = spiceypy.str2et(time_str)
    matrix = spiceypy.pxform(frame_in, frame_out, et)
    vt[i,:] = spiceypy.mxv(matrix, v[i,:])
  execution_stop = time.time()

  spiceypy.kclear()

  return vt, execution_stop - execution_start


def _spacepy(v, t, frame_in, frame_out, lib):
  import time
  import numpy as np

  import spacepy.coordinates as sc
  from spacepy.time import Ticktock

  import hxform

  import warnings
  warnings.filterwarnings("ignore", message="Leapseconds.*")

  if len(t.shape) == 1:
    # SpacePy requires time values to be strings with 1-second precision
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.ints_pad(t, length=6))
  else:
    t_str = []
    for i in range(t.shape[0]):
      t_str.append('%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.ints_pad(t[i,:], length=6)))
    t_str = np.array(t_str)

  execution_start = time.time()
  if lib.endswith('-irbem'):
    cvals = sc.Coords(v, frame_in, 'car', use_irbem=True)
  else:
    cvals = sc.Coords(v, frame_in, 'car', use_irbem=False)

  cvals.ticks = Ticktock(t_str, 'ISO')
  vt = cvals.convert(frame_out, 'car').data
  execution_stop = time.time()

  return vt, execution_stop - execution_start


def _sunpy(v, t, frame_in, frame_out):
  import time
  import hxform
  import numpy as np
  import astropy.coordinates
  import sunpy.coordinates
  # sunpy.coordinates is not used directly, but needed to register the frames.

  vt = np.full(v.shape, np.nan)

  representation_type = 'cartesian'

  # Does not seem possible to transform a dimensionless vector.
  one = astropy.constants.R_earth
  #one = astropy.units.m
  # TODO: Use
  #   one = astropy.units.m
  # when https://github.com/sunpy/sunpy/pull/7530 is merged.

  info = hxform.lib_info('sunpy')
  frame_in = info['system_aliases'][frame_in]
  frame_out = info['system_aliases'][frame_out]

  execution_start = time.time()
  if t.shape[0] == 1:
    obstime = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.ints_pad(t[0, :], length=6))
    kwargs = {
      "x": v[:, 0]*one,
      "y": v[:, 1]*one,
      "z": v[:, 2]*one,
      "frame": frame_in,
      "obstime": obstime,
      "representation_type": representation_type
    }
    coord = astropy.coordinates.SkyCoord(**kwargs)
    coord = coord.transform_to(frame_out).cartesian/one
    vt = coord.xyz.decompose().value.transpose()
  else:
    for i in range(t.shape[0]):
      obstime = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(hxform.timelib.ints_pad(t[i, :], length=6))
      kwargs = {
        "x": v[i, 0]*one,
        "y": v[i, 1]*one,
        "z": v[i, 2]*one,
        "frame": frame_in,
        "obstime": obstime,
        "representation_type": representation_type
      }

      coord_in = astropy.coordinates.SkyCoord(**kwargs)
      coord_out = coord_in.transform_to(frame_out).cartesian/one

      vt[i, :] = coord_out.xyz.decompose().value

  execution_stop = time.time()

  return vt, execution_stop - execution_start


def _sscweb(v, t, frame_in, frame_out):
  import re
  import requests
  import numpy as np
  import hxform
  import time as time
  vt = np.full(v.shape, np.nan)

  info = hxform.lib_info('sscweb')
  frame_in = info['system_aliases'].get(frame_in, frame_in)
  frame_out = info['system_aliases'].get(frame_out, frame_out)

  if frame_in in ['GEO', 'GM']:
    v[:,0], v[:,1], v[:,2] = car2sph(v[:,0], v[:,1], v[:,2])

  execution_start = time.time()
  for i in range(t.shape[0]):

    time.sleep(0.1) # To avoid overwhelming the server with requests and getting blocked.
    year = t[i][0]
    doy_ = hxform.timelib.doy(t[i][0:3])
    time_str = f'{year} {doy_:03d} {t[i][3]:02d}:{t[i][4]:02d}:{t[i][5]:02d}'
    # TODO: Document why SSCWeb Python client was not used
    #       (it handles too many request and retries)
    url = "https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi?"
    if frame_in in ['GEO', 'GM']:
      url += f"epoch={time_str}&x=&y=&z=&lat={v[i][1]:f}&lon={v[i][2]:f}&r={v[i][0]:f}&action={frame_in}"
    else:
      url += f"epoch={time_str}&x={v[i][0]:f}&y={v[i][1]:f}&z={v[i][2]:f}&lat=&lon=&r=&action={frame_in}"

    # This is likely to fail if many requests are made.
    try:
      print(f"Fetching URL: {url}")
      r = requests.get(url, timeout=2)
    except Exception as e:
      #print(e)
      raise Exception(f"Failed to fetch URL: {url}.")
      vt[i] = np.full((1, 3), np.nan)
      continue

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
      #logger.error(f"Failed to find start of table in URL: {url}. Returned HTML:\n{text}")
      vt[i] = np.full((1, 3), np.nan)
      continue

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

    vt[i][0] = result[frame_out]['X']
    vt[i][1] = result[frame_out]['Y']
    vt[i][2] = result[frame_out]['Z']

  execution_stop = time.time()

  return vt, execution_stop - execution_start


def matrix(t, frame_in, frame_out, lib='geopack_08_dp'):
  import numpy as np
  kwargs = {'ctype_in': 'car', 'ctype_out': 'car', 'lib': lib}
  c1 = transform(np.array([1., 0., 0.]), t, frame_in, frame_out, **kwargs)
  c2 = transform(np.array([0., 1., 0.]), t, frame_in, frame_out, **kwargs)
  c3 = transform(np.array([0., 0., 1.]), t, frame_in, frame_out, **kwargs)
  if len(c1.shape) == 1:
    m = np.column_stack([c1, c2, c3])
    if isinstance(t, list):
      return m.tolist()
    return m

  m = np.full((t.shape[0], 3, 3), np.nan)
  for i in range(t.shape[0]):
    m[i, :, 0] = c1[i, :]
    m[i, :, 1] = c2[i, :]
    m[i, :, 2] = c3[i, :]
  if isinstance(t, list):
    return m.tolist()
  return m


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