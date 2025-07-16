import datetime
import itertools

import numpy
import pandas

import hxform
import utilrsw

log = utilrsw.logger(console_format='%(message)s')

def libs(csys_in, csys_out, excludes=None):
  lib_infos = hxform.info.known_libs(info=True)
  libs_avail = []
  for i, lib_info in enumerate(lib_infos):
    lib = lib_info['name']
    if excludes is not None and lib in excludes:
      continue
    if csys_in not in lib_info['systems'] or csys_out not in lib_info['systems']:
      log.warning(f"Skipping {lib} because it does not support both {csys_in} and {csys_out} systems")
      continue
    libs_avail.append(lib)
  return libs_avail

def angles(to, tf, delta, libs, transform_kwargs):
  t = hxform.time_array(to, tf, delta)
  t_dts = hxform.ints2datetime(t)
  angles = numpy.full((t.shape[0], len(libs)), numpy.nan)
  for i, lib in enumerate(libs):
    #log.info(f"Processing {lib}...")
    transform_kwargs['lib'] = lib

    p_in = numpy.array([0., 0., 1.])
    p_out = hxform.transform(p_in, t, **transform_kwargs)

    n = numpy.dot(p_out, p_in)
    d = numpy.linalg.norm(p_out, axis=1)*numpy.linalg.norm(p_in)
    angles[:, i] = (180.0/numpy.pi)*numpy.arccos(n/d)

    if lib.startswith('spiceypy1'):
      years = numpy.array([dt.year for dt in t_dts])
      mask = (years < 1990) | (years >= 2020)
      angles[mask, i] = numpy.nan

    if lib.startswith('spiceypy2'):
      years = numpy.array([dt.year for dt in t_dts])
      mask = (years < 1990) | (years >= 2030)
      angles[mask, i] = numpy.nan

  return t_dts, angles


delta = {'days': 1}
to = datetime.datetime(2010, 1, 1, 0, 0, 0)
tf = datetime.datetime(2015, 1, 1, 0, 0, 0)
excludes = ['sscweb', 'cxform']

delta = {'days': 10}
to = datetime.datetime(2010, 1, 1, 0, 0, 0)
tf = datetime.datetime(2010, 12, 31, 0, 0, 0)
excludes = ['cxform']

delta = {'days': 1}
to = datetime.datetime(2010, 12, 21, 0, 0, 0)
tf = datetime.datetime(2010, 12, 23, 0, 0, 0)
excludes = ['cxform']

to_str = datetime.datetime.strftime(to, '%Y%m%d')
tf_str = datetime.datetime.strftime(tf, '%Y%m%d')
delta_unit = list(delta.keys())[0]
delta_str = f'{delta[delta_unit]}{delta_unit}'
file_out = f'data/dipole/dipole_delta={delta_str}_{to_str}-{tf_str}.pkl'

transform_kwargs = {
  'ctype_in': 'car',
  'ctype_out': 'car'
}

csys_list = ['GEO', 'MAG', 'GSE', 'GSM']
combinations = list(itertools.combinations(csys_list, 2))

dfs = {}
df_deltas = {}
for combination in combinations:
  csys_in, csys_out = combination
  #log.info(f"Processing {csys_in} to {csys_out}...")
  transform_kwargs['csys_in'] = csys_in
  transform_kwargs['csys_out'] = csys_out
  libs_avail = libs(csys_in, csys_out, excludes=excludes)
  t_dts, angles_ = angles(to, tf, delta, libs_avail, transform_kwargs)

  xform = f"{csys_in}_{csys_out}"
  dfs[xform] = {'values': None, 'diffs': None}
  dfs[xform]['values'] = pandas.DataFrame(angles_, columns=libs_avail)
  dfs[xform]['values'].index = pandas.to_datetime(t_dts)

  columns_diff = []
  for column in dfs[xform]['values'].columns:
    if column == 'geopack_08_dp':
      continue
    columns_diff.append(f'{column} diff')

  dfs[xform]['diffs'] = pandas.DataFrame(numpy.nan, index=dfs[xform]['values'].index, columns=columns_diff)
  if 'geopack_08_dp' in dfs[xform]['values'].columns:
    for lib in libs_avail:
      if lib != 'geopack_08_dp':
        col_name = f"{lib} diff"
        dfs[xform]['diffs'][col_name] = dfs[xform]['values'][lib] - dfs[xform]['values']['geopack_08_dp']

  dfs[xform]['diffs']['|max-min|'] = numpy.abs(dfs[xform]['values'].max(axis=1) - dfs[xform]['values'].min(axis=1))

  df_str = dfs[xform]['values'].to_string()
  log.info(f"{xform}\n{df_str}")
  df_str = dfs[xform]['diffs'].to_string()
  log.info(f"{xform}\n{df_str}")

print(f"Writing {file_out}")
utilrsw.write(file_out, dfs)

utilrsw.rm_if_empty('dipole.error.log')
