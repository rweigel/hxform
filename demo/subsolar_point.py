import datetime
import pandas

import hxform

import utilrsw

def compute(opts):

  t = utilrsw.time.ints_list(opts['to'], opts['tf'], opts['delta'])

  times = []
  thetas = []
  phis = []
  frame = opts['frame']
  for ti in t:
    times.append(datetime.datetime(*ti))

    theta_row = []
    phi_row = []
    for lib in opts['libs']:

      try:
        hxform.xprint(f'{lib} at {ti}')
        xyz = hxform.subsolar_point(ti, frame=opts['frame'], lib=lib)
        rtp = hxform.car2sph(xyz) # radius, theta (colat), phi (lon)
        theta_row.append(rtp[1])
        phi_row.append(rtp[2])

        hxform.xprint(f"  {frame} x, y, z: {xyz[0]:.6f}, {xyz[1]:.6f}, {xyz[2]:.6f}")
        hxform.xprint(f"  {frame} Lat, Lon: {rtp[1]:.6f}, {rtp[2]:.6f}")
      except Exception as e:
        theta_row.append(None)
        phi_row.append(None)

        hxform.xprint(f'  Error: {e}')

    thetas.append(theta_row)
    phis.append(phi_row)

  thetas = pandas.DataFrame(thetas, columns=opts['libs'], index=times)
  phis = pandas.DataFrame(phis, columns=opts['libs'], index=times)

  return thetas, phis


opts = {
  'to': [1965, 1, 1, 0, 0, 0],
  'tf': [2024, 1, 1, 0, 0, 0],
  'delta': {'hours': 24},
#  'delta': {'days': 50},
  'frame': 'GEO',
#  'libs': ['geopack_08_dp', 'cxform', 'spacepy'],
  'libs': None,
  'lib_ref': 'geopack_08_dp',
  'libs_exclude': ['sscweb']
}

# If no libraries specified, use all available
if opts['libs'] is None:
  libs_alt = ['tiegcm', 'apex', 'laundal', 'mead']
  opts['libs'] = [*libs_alt, *hxform.libs()]

# Remove excluded libraries
if opts['libs_exclude'] is not None:
  opts['libs'] = [lib for lib in opts['libs'] if lib not in opts['libs_exclude']]

thetas, phis = compute(opts)

data = {
  'opts': opts,
  'thetas': thetas,
  'phis': phis
}

with open('subsolar_point/subsolar_point.pkl', 'wb') as f:
  import pickle
  hxform.xprint("Writing subsolar_point/subsolar_point.pkl")
  pickle.dump(data, f)
