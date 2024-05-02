import datetime

import numpy as np
from hxform import hxform

import matplotlib.pyplot as plt

from datetick import datetick

def datetime2ints(dt):
  return [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]

if True:
  # 12:06:40+00:00 => 6/(1440) + 40/86400 = 0.00462962962
  # The time of the vernal equinox should be noon J2000.
  # Time: 2000-03-20T12:06:40.000000000
  # jd            utc   2451624.0046296297
  #                           0.00462962962 # From above
  # Vernal equinox is 2000-03-20T07:25 UTC
  n_steps = 1440
  delta_secs = 60
  geo = np.zeros((n_steps, 3))

  dt = datetime.datetime(1999, 12, 31, 12, 0, 0, tzinfo=datetime.timezone.utc)
  dtsi = []
  dts = []
  for i in range(n_steps):
    dtsi.append(datetime2ints(dt))
    dts.append(dt)
    dt = dt + datetime.timedelta(seconds=delta_secs)

  geo = hxform.transform(np.array([1, 0, 0]), dtsi, 'GEI', 'GEO', lib='spiceypy1')

  _, axes = plt.subplots(1, figsize=(8,4) )
  axes.plot(dts, geo[:,0])
  axes.grid()
  datetick()