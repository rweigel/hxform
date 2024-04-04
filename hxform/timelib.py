def tpad(time, length=7):
  """Pad list or array with 3 or more elements with zeros.

  Example:
  --------
  >>> from hxform import hxform as hx
  >>> print(hx.tpad([2000,1,1]))                 # [2000, 1, 1, 0, 0, 0, 0]
  >>> print(hx.tpad([2000,1,1], length=4))       # [2000, 1, 1, 0]
  >>> print(hx.tpad([2000,1,1,2,3,4], length=3)) # [2000, 1, 1]
  """
  import numpy as np
  in_type = type(time)

  # TODO: Check that time is valid
  time = np.array(time)

  assert len(time) > 2, "time must have at least 3 elements"

  if len(time.shape) == 1:
    if len(time) > length:
      time = time[0:length]
    else:
      pad = length - len(time)
      time = np.pad(time, (0, pad), 'constant', constant_values=0)
  else:
    if len(time[0]) > length:
      time = time[:,0:length]
    else:
      pad = length - len(time)
      time = np.pad(time, ((0, 0), (0, pad)), 'constant', constant_values=0)

  if in_type == np.ndarray:
    return time
  elif in_type == tuple:
    return tuple(map(tuple, time))
  else:
    return list(time)


def is_leap_year(year):
  import numpy as np
  if isinstance(year, int):
    if year % 100 == 0:
      return year % 400 == 0
    return year % 4 == 0
  else:
    year = np.array(year)
    leap400 = np.where(year % 400 == 0, True,False)
    leap100 = np.where(year % 100 == 0, True, False)
    leap4 = np.where(year %4 == 0, True, False)
    cor1 = np.where(leap100, False, leap4)
    return np.where(leap400, True, cor1)


def doy(date):
  """
  doy([2021,7,5])  ->  186
  doy([[2000,9,30],[1900,9,30],[1980,5,5]]) -> [274, 273. 126]
  """
  import numpy as np
  # https://astronomy.stackexchange.com/questions/2407/calculate-day-of-the-year-for-a-given-date
  date = np.array(date)

  if len(date.shape) == 1:
    year, month, day = date[0], date[1], date[2]

    if is_leap_year(year):
      K = 1
    else:
      K = 2
  else:
    year, month, day = date[:,0], date[:,1], date[:,2]
    K = np.where(is_leap_year(year),1, 2)

  N = np.fix((275.0*month)/9.0) - K*np.fix((month + 9.0)/12.0) + day - 30.0

  return N.astype(int)


def ints2doy(t):
  """Convert from [y, m, d, ...] to [y, doy, ...].

  Examples
  --------
  >>> ints2doy([2000, 2, 1, 9, 9, 9]) # [2000, 32, 9, 9, 9]
  """
  from datetime import datetime
  import numpy as np

  in_type = type(t)

  t = np.array(t)

  if len(t.shape) == 1:
    pad = 6 - len(t)
    t = np.pad(t, (0, pad), 'constant', constant_values=0)
  else:
    pad = 6 - len(t[0])
    t = np.pad(t, ((0, 0), (0, pad)), 'constant', constant_values=0)

  if len(t.shape) == 1:
    day_of_year = datetime(*t).timetuple().tm_yday
  else:
    day_of_year = doy(t[:,:3])
    t = np.column_stack((t[:,0], day_of_year, t[:,3], t[:,4], t[:,5]))

  if in_type == np.ndarray:
    return t
  elif in_type == tuple:
    return tuple(t.tolist()) #tuple(map(tuple,t))
  else:
    return t.tolist()


def iso2ints(isostr, length=None):
  """Convert time string in for YYYY-MM-DD[THH:mm:SS.FZ] to list of integers."""
  import re

  if not isinstance(isostr, str):
    if not isinstance(isostr, list):
      raise ValueError('Input must be a string or list of strings.')
    for i in range(len(isostr)):
      isostr[i] = iso2ints(isostr[i], length=7)
    return isostr

  tmp = re.split("-|:|T|Z", isostr)
  if len(tmp) > 6:
    tmp = tmp[0:5]

  int_list = []
  for str_int in tmp:
    if str_int != "Z" and str_int != '':
      int_list.append(int(str_int))

  return int_list


def UTtoHMS(UT, **kwargs):
  """Convert universal time in fractional hours into integer hour, minutes, seconds.

  Example
  -------
  >>> from hxform import hxform as hx
  >>> print(hx.UTtoHMS(12))              # [12, 0, 0]
  >>> print(hx.UTtoHMS(24))              # [0, 0, 0]
  >>> print(hx.UTtoHMS(24, keep24=True)) # [24, 0, 0]
  """

  keep24 = False
  if 'keep24' in kwargs:
    keep24 = kwargs['keep24']

  if UT > 24 or UT < 0:
    raise ValueError('Required: 0 <= UT <= 24.')

  hours = int(UT)
  minutes = int((UT-hours)*60.)
  seconds = int(round((UT-hours-minutes/60.)*3600.))
  if seconds == 60:
    seconds = 0
    minutes = minutes + 1
  if minutes == 60:
    minutes = 0
    hours = hours + 1

  if hours == 24 and keep24 == False:
    return [0, 0, 0]

  return [hours, minutes, seconds]
