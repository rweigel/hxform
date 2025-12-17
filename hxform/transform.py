def transform(v, t, frame_in, frame_out, ctype_in='car', ctype_out='car', lib='cxform'):
  """Transform between coordinates frames using cxform, Geopack-08, PySPEDAS, 
  SpacePy, SpiceyPy, SSCWeb, or SunPy.

  Parameters
  ----------
  v : array-like

      (Nv, 3) ndarrays, where `Nv` is the number of vectors to transform and
      v has dytpe of np.number. This option will result in the smallest execution time.

      list/tuple of three numbers (ints, floats, or numpy.numbers)

      list/tuple containing `Nv` lists/tuples of three numbers

      list/tuple of `Nv` (3, ) or (3, 1) ndarrays

  t : str or array-like

      In the following, `Nt` is the number of times to transform at. All times
      are interpreted as in UTC.

      str in the exact ISO format 'YYYY-MM-DDTHH:MM:SSZ' or 'YYYY-MM-DDTHH:MM:SS'.

      list/tuple of `Nt` strs in ISO format above and `Nt` = 1 or `Nt` = `Nv`.

      list/tuple of 3+ ints ([year, month, day, [hours, [minutes, [seconds]]]])
      Zeros are used for any missing optional value.

      list/tuple containing `Nt` lists/tuples of 3+ ints, where `Nt` = 1 or `Nt` = `Nv`.

      (Nt, 3+) int ndarray, where `Nt` = 1 or `Nt` = `Nv`.


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
  array-like with dimensions matching either `t` (if `Nt` ≠ 1 and `Nv` = 1) or
  `v` (if `Nv`≠ 1 and `Nt` = 1). If `Nv` and `Nt`≠ 1, dimensions are same as `v`.

  Use `np.ndarrays` for `v` and `t` and `t` with ints for fastest execution time when
  cxform and Geopack-08 libraries are used.

  Return type will match that of `v`. Note that if a list/tuple of 3-element np.arrays are
  passed, execution time will be long. 

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

  _arg_check(frame_in, frame_out, ctype_in, ctype_out, lib)

  # Check t and convert to standard form of (Nt, 6) int ndarray
  t, Nt = _t_prep(t)

  # Check v and convert to standard form of (Nv, 3) float ndarray
  v, v_outertype, v_innertype = _v_prep(v)

  Nv = v.shape[0]  # Number of vectors

  if Nt != 1 and Nv != 1 and Nt != Nv:
    msg = "If number of vectors and number of times are both > 1, they must be equal."
    raise ValueError(msg)

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


def _arg_check(frame_in, frame_out, ctype_in, ctype_out, lib):
  import hxform
  assert lib in hxform.libs(), f'lib must be one of {hxform.info.libs()}'
  assert ctype_in in ['car', 'sph'], 'ctype_in must be one of ["car", "sph"]'
  assert ctype_out in ['car', 'sph'], 'ctype_out must be one of ["car", "sph"]'

  info = hxform.lib_info(lib)
  for csys in [frame_in, frame_out]:
    emsg = f'For lib={lib}, {csys} must be one of {info["frames"]}'
    assert frame_in in info['frames'], emsg


def _v_prep(v):
  import numpy as np

  if not isinstance(v, (list, tuple, np.ndarray)):
    raise ValueError(f"v must be a list, tuple, or np.ndarray, not {type(v).__name__}")

  if isinstance(v, (list, tuple)) and len(v) == 0:
    raise ValueError(f"v is an empty {type(v).__name__}")

  if isinstance(v, np.ndarray):
    if v.ndim == 0:
      raise ValueError("v may not be a scalar (0-dim ndarray)")
    if v.size == 0:
      raise ValueError("v.size = 0")

  v_outertype = type(v)
  v_innertype = type(v[0])

  from utilrsw.np import components2matrix
  try:
    # Check v and convert to standard form of (Nv, 3) float ndarray
    v = components2matrix(v)
  except Exception as e:
    raise ValueError(f"utilrsw.np.components2matrix(v) raised '{e}'") from e

  return v, v_outertype, v_innertype


def _t_prep(t):

  import numpy as np

  def _t2ints(t):
    import datetime
    fmt = "%Y-%m-%dT%H:%M:%S"
    t_parsed = []
    for ti in t:
      if not isinstance(ti, str):
        raise ValueError("If time is a list, tuple, or ndarray of strings, all elements must be strings.")
      if ti.endswith('Z'):
        ti = ti[:-1]
      try:
        dt = datetime.datetime.strptime(ti, fmt)
      except ValueError as e:
        emsg = f"datetime.datetime.strptime('{ti}', '{fmt}') failed"
        raise ValueError(emsg) from e

      t_parsed.append([dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second])

    t = t_parsed
    if len(t) == 1:
      t = t[0]

    return t

  if isinstance(t, str):
    t = [t]

  if not isinstance(t, (list, tuple, np.ndarray)):
    raise ValueError("time must be a str, list, tuple, or np.ndarray")

  if isinstance(t, (list, tuple)) and len(t) == 0:
    raise ValueError("time input is empty")

  if isinstance(t, np.ndarray) and t.size == 0:
    raise ValueError("time input is empty")

  if isinstance(t, np.ndarray) and t.ndim == 0:
    t = np.array([t.item()])

  if isinstance(t[0], str):
    t = _t2ints(t)

  try:
    t = np.array(t, dtype=np.int32)
  except:
    raise ValueError("Invalid time input.")

  if len(t.shape) > 2:
    raise ValueError("Invalid time input.")

  if len(t.shape) == 1:
    t = np.array([t])

  prefix = "When represented as a list or tuple of ints, each time"
  if t.shape[1] > 6:
    raise ValueError(f"{prefix} must be given by 3 to 6 ints.")

  if t.shape[1] < 3:
    raise ValueError(f"{prefix} must have at least 3 elements.")

  if t.shape[1] < 6:
    n_pad = 6 - t.shape
    t = np.pad(t, ((0, 0), (0, n_pad)), 'constant', constant_values=0)

  Nt = t.shape[0]  # Number of times

  return t, Nt


def _cxform(v, t, frame_in, frame_out):
  import os
  import time
  import glob
  import ctypes

  import numpy as np

  lib_dir = os.path.join(os.path.dirname(__file__), "..", "hxform", "cxform_wrapper*")

  for lib_file in glob.glob(lib_dir):
    # The name of the .so or .dll file will not be the same on all
    # frames, so we need to find it. (For example, on one system
    # it is cxform_wrapper.cpython-37m-darwin.so.)
    # TODO: Find a better way to do this.
    break

  lib_path = os.path.join(lib_file)
  lib_obj = ctypes.cdll.LoadLibrary(lib_path)

  execution_start = time.time()

  # Allocate output array
  Nt = t.shape[0]
  if Nt == 1:
    vt = np.full(v.shape, np.nan)
  else:
    vt = np.full((Nt, 3), np.nan)

  # TODO: Handle ret
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

  import numpy

  import hxform
  from hxform import geopack_08_dp

  trans = frame_in + 'to' + frame_out
  dtime = numpy.array(hxform.timelib.ints2doy(t))

  # Allocate output array
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

  import numpy

  os.environ["PYSPEDAS_LOGGING_LEVEL"] = "error"
  from pyspedas import cotrans
  from pyspedas import time_double

  time_in = []
  for i in range(t.shape[0]):
    tstr = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(t[i, :])
    time_in.append(time_double(tstr))

  execution_start = time.time()
  vt = cotrans(time_in=time_in, data_in=v, coord_in=frame_in, coord_out=frame_out)
  if isinstance(vt, list):
    # https://github.com/spedas/pyspedas/issues/1273
    vt = numpy.array(vt)
  execution_stop = time.time()

  return vt, execution_stop - execution_start


def _spiceypy(v, t, frame_in, frame_out, lib):

  import os
  import time

  import numpy
  import spiceypy

  import hxform

  vt = numpy.full(v.shape, numpy.nan)

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
    time_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(t[i,:])
    et = spiceypy.str2et(time_str)
    matrix = spiceypy.pxform(frame_in, frame_out, et)
    vt[i,:] = spiceypy.mxv(matrix, v[i,:])
  execution_stop = time.time()

  spiceypy.kclear()

  return vt, execution_stop - execution_start


def _spacepy(v, t, frame_in, frame_out, lib):
  import time

  import numpy

  import spacepy.coordinates as sc
  from spacepy.time import Ticktock

  import warnings
  warnings.filterwarnings("ignore", message="Leapseconds.*")

  if len(t.shape) == 1:
    # SpacePy requires time values to be strings with 1-second precision
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(t)
  else:
    t_str = []
    for i in range(t.shape[0]):
      t_str.append('%04d-%02d-%02dT%02d:%02d:%02d' % tuple(t[i,:]))
    t_str = numpy.array(t_str)

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

  import numpy

  import astropy.coordinates
  import sunpy.coordinates # Not used directly, but needed to register the frames.

  vt = numpy.full(v.shape, numpy.nan)

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
    obstime = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(t[0, :])
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
      obstime = '%04d-%02d-%02dT%02d:%02d:%02d' % tuple(t[i, :])
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
  import time as time

  import numpy

  import hxform
  vt = numpy.full(v.shape, numpy.nan)

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
      vt[i] = numpy.full((1, 3), numpy.nan)
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
      vt[i] = numpy.full((1, 3), numpy.nan)
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
