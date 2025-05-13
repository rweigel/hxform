import os

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

import hapiplot
from hapiclient import hapi
from hapiclient import hapitime2datetime

# Earth radius in km used by SSCWeb. See
# https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
R_E = 6378.16 # km

################################################################################
# Options
################################################################################
#satellites = '_all_' 
#satellites = ['cluster', 'dscovr', 'geotail', 'themis', 'mms']
satellites = ['cluster']

#import sunpy
#sunpy.log.setLevel('DEBUG')

hapi_logging = True
interp_times = True
print_vals = False
print_first_last = True

legend_kwargs = {
  'loc': 'upper center',
  'markerscale': 3,
  'borderaxespad': 0,
  'framealpha': 1,
  'frameon': True,
  'ncol': 3
}

matplotlib.use('Agg')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 14
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.figsize'] = (8.5, 11)
################################################################################

def infos(satellite, source):

  """
  SSCWeb provides TOD, J2K, GEO, GM, GSE, and SM
    TOD in SSCWeb HAPI server is same as GEI from SSCWeb.
    (The reason the HAPI server uses TOD for GEI is that "TOD" is the 
    POST query parameter names used to request data in GEI is "TOD" and not GEI
    and consistency in the query parameters was sought).
  """

  infos_cdaweb = []
  infos_sscweb = []
  infos_amda = []

  if satellite == '' or satellite == 'dscovr':
    start = '2021-11-25T00:00:00Z'
    stop  = '2021-12-05T00:00:00Z'

    for frame in ['GSE', 'GCI']:
      infos_cdaweb.append({
        'dataset': 'DSCOVR_ORBIT_PRE',
        'parameter': f'{frame}_POS',
        'start': start,
        'stop':  stop
      })

      if frame == 'GCI':
        frame = 'J2K'

      infos_sscweb.append({
          'dataset': 'dscovr',
          'frame': frame,
          'start': start,
          'stop':  stop
        })

    infos_cdaweb.append(infos_cdaweb[-1].copy())
    infos_sscweb.append(infos_sscweb[-1].copy())
    infos_sscweb[-1]['frame'] = 'GEI'

  if satellite == '' or satellite == 'themis':

    themis = 'themisa'
    start = '2016-09-01T00:00:00Z'
    stop  = '2016-09-02T00:00:00Z'
    for frame in ['GSE', 'GSM', 'GEI']:
      if frame == 'GEI':
        parameter = f'th{themis[-1]}_pos'
      else:
        parameter = f'th{themis[-1]}_pos_{frame.lower()}'

      infos_cdaweb.append({
        'dataset': f'TH{themis[-1].upper()}_L1_STATE@0',
        'parameter': parameter,
        'start': start,
        'stop':  stop
      })

      if frame == 'GEI':
        frame = 'J2K'

      infos_sscweb.append({
        'dataset': themis,
        'frame': frame,
        'start': start,
        'stop':  stop
      })

    infos_cdaweb.append(infos_cdaweb[-1].copy())
    infos_sscweb.append(infos_sscweb[-1].copy())
    infos_sscweb[-1]['frame'] = 'TOD'

  if satellite == '' or satellite == 'mms':

    mms   = 'mms1'
    start = '2016-09-01T00:00:00Z'
    stop  = '2016-09-02T00:00:00Z'
    for frame in ['GSE', 'GSM']:
      infos_cdaweb.append({
        'dataset': f'{mms.upper()}_EPD-EIS_SRVY_L2_ELECTRONENERGY',
        'parameter': f'mms1_epd_eis_srvy_l2_electronenergy_position_{frame.lower()}',
        'start': start,
        'stop':  stop
      })

      infos_sscweb.append({
        'dataset': mms,
        'frame': frame,
        'start': start,
        'stop':  stop
      })

  if satellite == '' or satellite == 'geotail':
    """
    From https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html, 
    "Geocentric Equatorial Inertial system. This system has X-axis
    pointing from the Earth toward the first point of Aries (the position of
    the Sun at the vernal equinox). This direction is the intersection of the
    Earth's equatorial plane and the ecliptic plane and thus the X-axis lies
    in both planes. The Z-axis is parallel to the rotation axis of the Earth,
    and y completes the right-handed orthogonal set (Y = Z * X). Geocentric
    Inertial (GCI) and Earth-Centered Inertial (ECI) are the same as GEI."

    Based on the above quote, we treat GCI from CDAWeb as the same as GEI from SSCWeb.

    https://cdaweb.gsfc.nasa.gov/misc/NotesG.html#GE_OR_DEF
    CDAWeb dataset has GCI, GSE, GSM, HEC
    """

    start      = '2021-11-25T00:00:00Z'
    stop       = '2021-12-05T00:00:00Z'

    for frame in ['GSE', 'GSM', 'GCI']:
      infos_cdaweb.append({
        'dataset': 'GE_OR_DEF',
        'parameter': f'{frame}_POS',
        'start': start,
        'stop':  stop
      })

      if frame == 'GCI':
        frame = 'J2K'

      infos_sscweb.append({
        'dataset': 'geotail',
        'frame': frame,
        'start': start,
        'stop':  stop
      })

    infos_cdaweb.append(infos_cdaweb[-1].copy())
    infos_sscweb.append(infos_sscweb[-1].copy())
    infos_sscweb[-1]['frame'] = 'GEI'

  if satellite == '' or satellite == 'cluster':

    sc = '1'
    dataset = f'C{sc}_CP_FGM_5VPS'
    start = '2010-09-01T09:09:00.100Z'
    stop  = '2010-09-02T00:00:00.000Z'

    for frame in ['GSE']:
      infos_cdaweb.append({
        'dataset': dataset,
        'parameter': f'sc_pos_xyz_{frame.lower()}__{dataset}',
        'start': start,
        'stop':  stop
      })

      infos_sscweb.append({
        'dataset': f'cluster{sc}',
        'frame': frame,
        'start': start,
        'stop':  stop
      })

      infos_amda.append({
        'dataset': f'clust{sc}-orb-all',
        'parameter': f'c{sc}_xyz_{frame.lower()}',
        'start': start,
        'stop':  stop
      })

    if source == 'sscweb':
      return infos_sscweb
    if source == 'cdaweb':
      return infos_cdaweb
    if source == 'amda':
      return infos_amda

def jpl(info_sscweb, logging=False):
  import pickle

  from sunpy.coordinates import get_horizons_coord

  # No 'GSM' b/c https://github.com/sunpy/sunpy/issues/8188
  if info_sscweb['frame'] not in ['GSE', 'GEI']:
    return None

  # To find ids, see https://ssd.jpl.nasa.gov/horizons/app.html#/
  # and edit the "Target Body" field to find the id. Does not seem possible
  # to query for list of Target Body ids.
  # Note: No Geotail or themis{a,d}
  known_ids = {
    'ace': -92,
    'dscovr': -78,
    'mms1': -140482,
    'mms2': -140483,
    'mms3': -140484,
    'mms4': -140485,
    'themisb': -192,
    'themisc': -193
  }

  if info_sscweb['dataset'] not in known_ids:
    return None
  satellite_id = known_ids[info_sscweb['dataset']]

  to = info_sscweb['time'][0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  tf = info_sscweb['time'][-1].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  cache_file = f"jpl/jpl-{info_sscweb['dataset']}-{info_sscweb['frame']}-{to}-{tf}.pkl"

  if False and os.path.exists(cache_file):
    print(f"  Saved pickle file: {cache_file}")
    with open(cache_file, 'rb') as f:
      return pickle.load(f)

  #solar_system_ephemeris.set('de432s')
  #solar_system_ephemeris.set('de440s')
  xyz = np.full((len(info_sscweb['time']), 3), np.nan)
  dt = info_sscweb['time'][1] - info_sscweb['time'][0]
  # Using, e.g., 60s gives error: https://github.com/sunpy/sunpy/issues/8188
  step_min = int(dt.total_seconds()/60)
  t = {
    "start": info_sscweb['time'][0],
    "stop": info_sscweb['time'][-1],
    "step": f'{step_min}m' 
  }

  if info_sscweb['frame'] == 'GSE':
    data_jpl = get_horizons_coord(satellite_id, t).geocentricsolarecliptic.cartesian.xyz.to('km')
  if info_sscweb['frame'] == 'GEI':
    data_jpl = get_horizons_coord(satellite_id, t).geocentricearthequatorial.cartesian.xyz.to('km')
  if info_sscweb['frame'] == 'GSM':
    data_jpl = get_horizons_coord(satellite_id, t).geocentricsolarmagnetospheric.cartesian.xyz.to('km')

  xyz = data_jpl.value.T/R_E
  with open(cache_file, 'wb') as f:
    pickle.dump({'time': info_sscweb['time'], 'xyz': info_sscweb['xyz']}, f)
    print(f"  Saved pickle file: {cache_file}")

  return {'time': info_sscweb['time'], 'xyz': xyz}

def sscweb(info, logging=False):

  server     = 'http://hapi-server.org/servers/SSCWeb/hapi'
  dataset    = info['dataset']
  frame      = info['frame'].replace('GEI', 'TOD')
  start      = info['start']
  stop       = info['stop']
  parameters = f'X_{frame},Y_{frame},Z_{frame}'

  opts       = {'logging': logging, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameters, start, stop, **opts)

  info['xyz'] = np.column_stack((data[f'X_{frame}'], data[f'Y_{frame}'], data[f'Z_{frame}']))

  # Convert from YYYY-DOY to YYYY-MM-DD date format
  info['time'] = hapitime2datetime(data['Time'])

  # Return not needed b/c info is modified in place. Keep for clarity.
  return info

def cdaweb(info, logging=False):

  server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
  dataset    = info['dataset']
  parameter  = info['parameter']
  start      = info['start']
  stop       = info['stop']
  opts       = {'logging': logging, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameter, start, stop, **opts)

  xyz = data[parameter]
  info['xyz'] = xyz/R_E # Convert from km to R_E
  # THe HAPI SSCWeb server provides ephemeris as three separate parameters. Here
  # we combine parameters into a list. SSCWeb reports in R_E while CDAWeb in km.
  # Convert CDAWeb data to R_E.
  info['time'] = hapitime2datetime(data['Time'])

  # Return not needed b/c info is modified in place. Keep for clarity.
  return info

def amda(info, logging=False):

  server     = 'https://amda.irap.omp.eu/service/hapi'
  dataset    = info['dataset']
  parameter  = info['parameter']
  start      = info['start']
  stop       = info['stop']
  opts       = {'logging': logging, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameter, start, stop, **opts)

  # https://amda.irap.omp.eu/service/hapi/info?id=clust1-orb-all
  # the referenced SPASE record says
  #"Units": "Re",
  #"UnitsConversion": "6.4e6>m"
  xyz = data[parameter]
  info['xyz'] = xyz*6.4e3/R_E # Convert from AMDA Re to SSCWeb R_E
  #info['xyz'] = xyz # Gives better match to CDAWeb
  info['time'] = hapitime2datetime(data['Time'])

  # Return not needed b/c info is modified in place. Keep for clarity.
  return info

def print_first_last_(info_cdaweb, info_sscweb, info_amda=None, info_jpl=None):
  if not print_first_last:
    return

  frame = info_sscweb['frame']
  for i in [0, -1]:
    xyz_cdaweb0 = info_cdaweb['xyz'][i]
    xyz_sscweb0 = info_sscweb['xyz'][i]
    if info_jpl is not None:
      xyz_jpl0 = info_jpl['xyz'][i]
    if info_amda is not None:
      xyz_amda0 = info_amda['xyz'][i]
    time_cdaweb0 = info_cdaweb['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    time_sscweb0 = info_sscweb['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if info_jpl is not None:
      time_jpl0 = info_jpl['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if info_amda is not None:
      time_amda0 = info_amda['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    u = 'R_E'
    if i == 0:
      print(f'  First values in {u}:')
    else:
      print(f'  Last values in {u}:')
    print('            {0:7s} {1:7s} {2:7s}'.format(f'X_{frame}', f'Y_{frame}', f'Z_{frame}'))
    print('    CDAWeb: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_cdaweb0, time_cdaweb0))
    print('    SSCWeb: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_sscweb0, time_sscweb0))
    if info_jpl is not None:
      print('       JPL: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_jpl0, time_jpl0))
    if info_amda is not None:
      print('      AMDA: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_amda0, time_amda0))

def interp(timei, time, xyz):
  # Interpolate CDAWeb data onto the union of timestamps
  interp_func = interp1d(
    [t.timestamp() for t in time],
    xyz,
    axis=0,
    bounds_error=False,
    fill_value="extrapolate"
  )

  return interp_func([t.timestamp() for t in timei])

def compute_diffs(info_cdaweb, info_sscweb):

  if interp_times:
    print("  Computing union of timestamps.")
    # Union of timestamps from CDAWeb and SSCWeb
    union_times = set(info_cdaweb['time']).union(set(info_sscweb['time']))
    union_times = sorted(union_times)

    print(f"  Interpolating to {len(union_times)} timestamps.")
    info_cdaweb['xyz'] = interp(union_times, info_cdaweb['time'], info_cdaweb['xyz'])
    info_cdaweb['time'] = union_times
    info_sscweb['xyz'] = interp(union_times, info_sscweb['time'], info_sscweb['xyz'])
    info_sscweb['time'] = union_times
  else:
    # Common timestamps
    print("  Computing common timestamps.")
    common_times = set(info_cdaweb['time']).intersection(set(info_sscweb['time']))
    common_times = sorted(common_times)
    common_times = sorted(common_times)
    print(f"  {len(common_times)} common timestamps.")

    # Find the indices of common times
    idx_cdaweb = np.nonzero(np.isin(info_cdaweb['time'], common_times))[0]
    idx_sscweb = np.nonzero(np.isin(info_sscweb['time'], common_times))[0]

    info_cdaweb['xyz'] = np.array(info_cdaweb['xyz'][idx_cdaweb])
    info_cdaweb['time'] = common_times
    info_sscweb['xyz'] = np.array(info_sscweb['xyz'][idx_sscweb])
    info_sscweb['time'] = common_times

  print("  Computing differences")
  #r_cdaweb = np.linalg.norm(info_cdaweb['xyz'], axis=1)
  #r_diff = np.diff(r_cdaweb)/R_E
  #t_diff = info_cdaweb['time'][1:]

  nt = len(info_cdaweb['time'])
  Δθ = np.full(nt, np.nan)
  Δr = np.full(nt, np.nan)
  Δr_rel = np.full(nt, np.nan)
  t = np.empty(nt, dtype='datetime64[ns]')
  t[:] = np.datetime64('NaT')
  r = compute_r(info_cdaweb, info_sscweb)

  info_cdaweb['time'] = [t.replace(tzinfo=None) for t in info_cdaweb['time']]

  for i in range(nt):
    # Remove timezone from info_cdaweb['time']
    t[i] = np.datetime64(info_cdaweb['time'][i])
    Δr[i] = np.linalg.norm(info_cdaweb['xyz'][i] - info_sscweb['xyz'][i])

    Δr_rel[i] = Δr[i]/r[i]

    # d = denominator for angle calculation
    # Δθ[i] = arccos [ (a·b)/(|a|*|b|) ] = arccos (n/d)
    n = np.dot(info_cdaweb['xyz'][i], info_sscweb['xyz'][i])
    d = np.linalg.norm(info_cdaweb['xyz'][i])*np.linalg.norm(info_sscweb['xyz'][i])
    if np.abs(n) > np.abs(d):
      if np.abs(n - d) > 1e-10:
        raise ValueError(f"n > d. n = {n:.16f} d = {d:.16f}, n-d = {n-d:.16f}")
      wmsg = f"  Warning: {t[i]}:\n    |n| > |d| in Δθ[i] = arccos [ (a·b)/(|a|*|b|) ]"
      wmsg += f" = arccos (n/d)\n    n = {n:.16f}, n-d = {n-d:.16f}"
      print(wmsg)
      # If n == d, then the angle is 0 degrees, but numpy gives
      # RuntimeWarning: invalid value encountered in arccos. Did rounding
      # cause n > d?
      Δθ[i] = 0
    else:
      Δθ[i] = (180/np.pi)*np.arccos(n/d)

    if print_vals:
      t_str = info_cdaweb['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
      print(f"t = {t_str} | Δr = {Δr[i]:.5f} [R_E] | Δ∠: {Δθ[i]:.5f}°")

  return t, Δr, Δr_rel, Δθ

def adjust_y_range(ax, gap_fraction=0.25, bottom=None):
  # Adjust the y-axis range to add a gap at the top so legend does not overlap
  # with the data. Seems that the needed gap could be computed.
  ylim = ax.get_ylim()
  yticks = ax.get_yticks()
  gap = gap_fraction * (yticks[1] - yticks[0])
  if bottom is None:
    bottom = ylim[0]
  ax.set_ylim(bottom, ylim[1] + gap)

def compute_r(info_cdaweb, info_sscweb):
  ra = np.linalg.norm(info_cdaweb['xyz'], axis=1)
  rb = np.linalg.norm(info_sscweb['xyz'], axis=1)
  # TODO: Interpolation to common times
  #return (ra + rb)/2
  return ra

def plot_xyz(ax, info_cdaweb, info_sscweb, info_jpl=None):
  colors = ['r', 'g', 'b']
  component = ['X', 'Y', 'Z']

  for c in range(3):
    ax.plot(info_cdaweb['time'], info_cdaweb['xyz'][:,c],
            label=f'CDAWeb/{component[c]}', lw=2, linestyle='-', color=colors[c])
    ax.plot(info_sscweb['time'], info_sscweb['xyz'][:,c],
            label=f'SSCWeb/{component[c]}', lw=3, linestyle='--', color=colors[c])
    if info_jpl is not None:
      ax.plot(info_jpl['time'], info_jpl['xyz'][:,c],
              label=f'JPL/{component[c]}', lw=1, linestyle=':', color=colors[c])

  r = compute_r(info_cdaweb, info_sscweb)
  ax.plot(info_cdaweb['time'], r,
          label='$\\overline{r}$', lw=2, linestyle='-', color='k')

  adjust_y_range(ax, gap_fraction=1)
  ax.set_ylabel('$R_E$', rotation=0)
  ax.grid()

  ax.legend(**{**legend_kwargs, 'ncol': 4})
  ax.set_xticklabels([])

def plot_diffs(ax, info_cdaweb, info_sscweb):
  print("  Computing differences")
  t, Δr, Δr_rel, Δθ  = compute_diffs(info_cdaweb, info_sscweb)

  lw = 2 # line width

  Δr_rel_max = np.nanmax(Δr_rel)
  #A = f'{Δr_rel_max:.1e}'.split('e')[0]
  #e = str(int(np.floor(np.log10(Δr_rel_max))))
  #Δr_rel_max_str = f" ($Δr_{{\\rm{{max}}}}/\\overline{{r}}$ = {A}·$10^{{{e}}}$)"
  Δr_rel_max_str = f" (max= 1/{1/Δr_rel_max:.0f})"

  Δr_max = np.nanmax(Δr)
  Δr_max_str = f'($Δr_{{\\rm{{max}}}} = {Δr_max*R_E:.1f}$ [km])'

  print(f"  Δr_max = {Δr_max:.5f} [R_E]")
  print(f"  Δr_max = {Δr_max*R_E:.1f} [km]")

  ax.plot(t, Δr, 'r-', lw=lw, label=f'$|Δr|/R_E$ {Δr_max_str}')
  ax.plot(t, Δr_rel, 'b-', lw=lw, label=f'$|Δr|/\\overline{{r}}$ {Δr_rel_max_str}')
  ax.plot(t, Δθ, 'g-', lw=lw, label='$Δθ$ [deg]')

  adjust_y_range(ax, bottom=0, gap_fraction=1)
  ax.legend(**legend_kwargs)
  ax.grid()

def _annotate(axes):
  axes[0].set_title(title)

  for ax in axes:
    ax.spines['bottom'].set_visible(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    # Remove minor tick lines next to y-axis numbers
    ax.tick_params(axis='y', which='minor', length=0)
    # Remove tick lines next to y-axis numbers
    ax.tick_params(axis='y', length=0)
    ax.grid(which='minor', linestyle=':', linewidth=0.5, color=[0.75]*3)
    ax.minorticks_on()

def _savefigs(dir, dataset, frame):
  for dir in ['svg', 'png', 'pdf']:
    os.makedirs(f'figures/{dir}', exist_ok=True)
    fname = f"figures/{dir}/{infos_sscweb[i]['dataset']}_{frame}"
    print(f"  Writing figures/{dir}/{fname}.{dir}")
    kwargs = {'bbox_inches': 'tight'}
    if dir == 'png':
      kwargs['dpi'] = 300
    plt.savefig(f'{fname}.{dir}', bbox_inches='tight')
  plt.close()

for satellite in satellites:

  infos_cdaweb = infos(satellite, 'cdaweb')
  infos_sscweb = infos(satellite, 'sscweb')
  infos_amda = infos(satellite, 'amda')
  print(infos_amda)
  for i in range(len(infos_cdaweb)):
    frame = infos_sscweb[i]['frame']

    a = f"SSCWeb/{infos_sscweb[i]['dataset']}/{frame}"
    b = f"CDAWeb/{infos_cdaweb[i]['dataset'].split('@')[0]}/{infos_cdaweb[i]['parameter']}"
    if infos_sscweb[i]['dataset'].startswith('mms'):
      # Shorten name for MMS
      b = f"CDAWeb/{infos_cdaweb[i]['parameter']}"
    print(f"Comparing\n  {a}\n  with\n  {b}")
    title = f"{a} vs. {b}"

    # Adds xyz and time to info dict
    infos_sscweb[i] = sscweb(infos_sscweb[i], logging=hapi_logging)
    infos_cdaweb[i] = cdaweb(infos_cdaweb[i], logging=hapi_logging)
    infos_amda[i] = amda(infos_amda[i], logging=hapi_logging)

    info_jpl = jpl(infos_sscweb[i], logging=hapi_logging)
    print_first_last_(infos_cdaweb[i], infos_sscweb[i], info_amda=infos_amda[i], info_jpl=info_jpl)

    gs = plt.gcf().add_gridspec(2)
    axes = gs.subplots(sharex=True)
    _annotate(axes)

    print("  Plotting xyz")
    plot_xyz(axes[0], infos_cdaweb[i], infos_sscweb[i], info_jpl=info_jpl)
    print("  Plotting diffs")
    plot_diffs(axes[1], infos_cdaweb[i], infos_sscweb[i])

    from datetime import timedelta
    # Subtract 1 minute to the x-axis limits to make sure the first tick is shown
    m = timedelta(minutes=1)
    # Round to the next hour + 1 minute
    hm = timedelta(hours=1, minutes=1)
    last = infos_cdaweb[i]['time'][-1].replace(minute=0, second=0, microsecond=0) + hm
    axes[0].set_xlim(infos_cdaweb[i]['time'][0] - m, last)
    axes[1].set_xlim(infos_cdaweb[i]['time'][0] - m, last)
    hapiplot.plot.datetick.datetick(dir='x')

    labels = axes[1].get_xticklabels()
    if labels:
      labels[0].set_horizontalalignment('left')
      labels[-1].set_horizontalalignment('right')

    _savefigs(dir, infos_sscweb[i]['dataset'], frame)

  print("")
