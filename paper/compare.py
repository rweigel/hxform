import os

import numpy as np
import matplotlib
from matplotlib import pyplot as plt

from hapiclient import hapi
from hapiclient import hapitime2datetime


"""
Earth radius in km used by SSCWeb. See
https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
"""
R_E = 6378.16 # km
"""
In https://amda.irap.omp.eu/service/hapi/info?id=clust1-orb-all
the referenced SPASE record says
"Units": "Re",
"UnitsConversion": "6.4e6>m"
"""
R_E_AMDA = 6.4e3 # km

################################################################################
# Options
################################################################################
satellites_only = [] # Run all satellites
#satellites_only = ['dscovr']

hapi_logging = False
cos_warnings = False
interp_times = True
print_first_last = True   # Print first and last uninterpolated values
print_interp_vals = False # Print interpolated values

import warnings
warnings.filterwarnings("ignore", message="The argument 'infer_datetime_format'")
warnings.filterwarnings("ignore", message="Could not infer format")

# Set to ERROR to suppress INFO message for making request to Horizons server
import sunpy
sunpy.log.setLevel('ERROR')

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

def infos(satellite=None):

  satellites = ['Cluster', 'DSCOVR', 'Geotail', 'THEMIS', 'MMS']
  if satellite is None:
    return satellites

  if satellite not in satellites:
    raise ValueError(f"Unknown satellite: {satellite}. Must be one of {satellites}")

  infos_ = {'cdaweb': [], 'sscweb': [], 'amda': []}

  if satellite == '' or satellite == 'Cluster':

    sc = '1'
    dataset = f'C{sc}_CP_FGM_5VPS'
    start = '2010-09-01T09:09:00.100Z'
    stop  = '2010-09-02T00:00:00.000Z'

    for frame in ['GSE', 'GSM']:

      if frame == 'GSE':
        infos_['cdaweb'].append({
          'name': 'CDAWeb',
          'dataset': dataset,
          'parameter': f'sc_pos_xyz_{frame.lower()}__{dataset}',
          'start': start,
          'stop':  stop
        })
      else:
        # No GSM data in CDAWeb
        infos_['cdaweb'].append(None)

      infos_['sscweb'].append({
        'name': 'SSCWeb',
        'dataset': f'cluster{sc}',
        'frame': frame,
        'start': start,
        'stop':  stop
      })

      infos_['amda'].append({
        'name': 'AMDA',
        'dataset': f'clust{sc}-orb-all',
        'parameter': f'c{sc}_xyz_{frame.lower()}',
        'start': start,
        'stop':  stop
      })

  if satellite == '' or satellite == 'DSCOVR':
    start = '2021-11-25T00:00:00Z'
    stop  = '2021-12-05T00:00:00Z'

    for frame in ['GSE', 'GCI']: # GCI must be last
      infos_["cdaweb"].append({
        'name': 'CDAWeb',
        'dataset': 'DSCOVR_ORBIT_PRE',
        'parameter': f'{frame}_POS',
        'start': start,
        'stop':  stop
      })

      infos_["sscweb"].append({
          'name': 'SSCWeb',
          'dataset': 'dscovr',
          'frame': frame,
          'start': start,
          'stop':  stop
        })

    infos_["cdaweb"].append(infos_["cdaweb"][-1].copy())
    infos_["sscweb"].append(infos_["sscweb"][-1].copy())
    infos_["sscweb"][-1]['frame'] = 'J2K'

  if satellite == '' or satellite == 'THEMIS':

    themis = 'themisa'
    start = '2016-09-01T00:00:00Z'
    stop  = '2016-09-02T00:00:00Z'
    for frame in ['GSE', 'GSM', 'GEI']: # GEI must be last
      if frame == 'GEI':
        # Metadata indicates COORDINATE_SYSTEM = GEI for this parameter
        parameter = f'th{themis[-1]}_pos'
      else:
        parameter = f'th{themis[-1]}_pos_{frame.lower()}'

      infos_['cdaweb'].append({
        'name': 'CDAWeb',
        'dataset': f'TH{themis[-1].upper()}_L1_STATE@0',
        'parameter': parameter,
        'start': start,
        'stop':  stop
      })

      infos_['sscweb'].append({
        'name': 'CDAWeb',
        'dataset': themis,
        'frame': frame,
        'start': start,
        'stop':  stop
      })

    infos_['cdaweb'].append(infos_['cdaweb'][-1].copy())
    infos_['sscweb'].append(infos_['sscweb'][-1].copy())
    infos_['sscweb'][-1]['frame'] = 'J2K'

  if satellite == '' or satellite == 'MMS':

    mms   = 'mms1'
    start = '2016-09-01T00:00:00Z'
    stop  = '2016-09-02T00:00:00Z'
    for frame in ['GSE', 'GSM']:
      infos_['cdaweb'].append({
        'name': 'CDAWeb',
        'dataset': f'{mms.upper()}_EPD-EIS_SRVY_L2_ELECTRONENERGY',
        'parameter': f'mms1_epd_eis_srvy_l2_electronenergy_position_{frame.lower()}',
        'start': start,
        'stop':  stop
      })

      infos_['sscweb'].append({
        'name': 'SSCWeb',
        'dataset': mms,
        'frame': frame,
        'start': start,
        'stop':  stop
      })

  if satellite == '' or satellite == 'Geotail':
    """
    https://cdaweb.gsfc.nasa.gov/misc/NotesG.html#GE_OR_DEF
    CDAWeb dataset has GCI, GSE, GSM, HEC
    """

    start      = '2021-11-25T00:00:00Z'
    stop       = '2021-12-05T00:00:00Z'

    for frame in ['GSE', 'GSM', 'GCI']: # GCI must be last
      infos_['cdaweb'].append({
        'name': 'CDAWeb',
        'dataset': 'GE_OR_DEF',
        'parameter': f'{frame}_POS',
        'start': start,
        'stop':  stop
      })

      infos_['sscweb'].append({
        'name': 'SSCWeb',
        'dataset': 'geotail',
        'frame': frame,
        'start': start,
        'stop':  stop
      })

    # Also compare GCI from CDAWeb with J2K from SSCWeb
    infos_['cdaweb'].append(infos_['cdaweb'][-1].copy())
    infos_['sscweb'].append(infos_['sscweb'][-1].copy())
    infos_['sscweb'][-1]['frame'] = 'J2K'

  if len(infos_['amda']) == 0:
    infos_['amda'] = len(infos_['cdaweb'])*[None]

  return infos_

def jpl(info, logging=False):
  import pickle

  from sunpy.coordinates import get_horizons_coord

  # No 'GSM' b/c https://github.com/sunpy/sunpy/issues/8188
  if info['frame'] not in ['GSE', 'GEI']:
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

  if info['dataset'] not in known_ids:
    return None
  satellite_id = known_ids[info['dataset']]

  to = info['time'][0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  tf = info['time'][-1].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  cache_file = f"jpl/jpl-{info['dataset']}-{info['frame']}-{to}-{tf}.pkl"

  if False and os.path.exists(cache_file):
    print(f"  Saved pickle file: {cache_file}")
    with open(cache_file, 'rb') as f:
      return pickle.load(f)

  #solar_system_ephemeris.set('de432s')
  #solar_system_ephemeris.set('de440s')
  xyz = np.full((len(info['time']), 3), np.nan)
  dt = info['time'][1] - info['time'][0]
  # Using, e.g., 60s gives error: https://github.com/sunpy/sunpy/issues/8188
  step_min = int(dt.total_seconds()/60)
  t = {
    "start": info['time'][0],
    "stop": info['time'][-1],
    "step": f'{step_min}m'
  }

  if info['frame'] == 'GSE':
    data_jpl = get_horizons_coord(satellite_id, t).geocentricsolarecliptic.cartesian.xyz.to('km')
  if info['frame'] == 'GEI':
    data_jpl = get_horizons_coord(satellite_id, t).geocentricearthequatorial.cartesian.xyz.to('km')
  if info['frame'] == 'GSM':
    data_jpl = get_horizons_coord(satellite_id, t).geocentricsolarmagnetospheric.cartesian.xyz.to('km')

  xyz = data_jpl.value.T/R_E
  with open(cache_file, 'wb') as f:
    pickle.dump({'time': info['time'], 'xyz': info['xyz']}, f)
    print(f"  Saved pickle file: {cache_file}")

  return {'name': 'JPL', 'time': info['time'], 'xyz': xyz}

def sscweb(info, logging=False):

  """
  1. SSCWeb provides TOD, J2K, GEO, GM, GSE, and SM
  TOD in SSCWeb HAPI server is same as GEI from SSCWeb.
  (The reason the HAPI server uses TOD for GEI is that "TOD" is the 
  POST query parameter names used to request data in GEI is "TOD" and not GEI
  and consistency in the query parameters was sought).

  2. From https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html,

  > "Geocentric Equatorial Inertial system. This system has X-axis
  pointing from the Earth toward the first point of Aries (the position of
  the Sun at the vernal equinox). This direction is the intersection of the
  Earth's equatorial plane and the ecliptic plane and thus the X-axis lies
  in both planes. The Z-axis is parallel to the rotation axis of the Earth,
  and y completes the right-handed orthogonal set (Y = Z * X). Geocentric
  Inertial (GCI) and Earth-Centered Inertial (ECI) are the same as GEI."

  Based on the above quote, we treat GCI from CDAWeb as the same as GEI from SSCWeb.
  """
  server     = 'http://hapi-server.org/servers/SSCWeb/hapi'
  dataset    = info['dataset']
  frame      = info['frame']
  if frame == 'GEI':
    frame = 'TOD' # See SSCWeb note 1. above.
  if frame == 'GCI':
    frame = 'TOD' # See SSCWeb note 2. above.
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

  if info is None:
    return None

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

  if info is None:
    return None

  server     = 'https://amda.irap.omp.eu/service/hapi'
  dataset    = info['dataset']
  parameter  = info['parameter']
  start      = info['start']
  stop       = info['stop']
  opts       = {'logging': logging, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameter, start, stop, **opts)

  xyz = data[parameter]
  info['xyz'] = xyz*R_E_AMDA/R_E # Convert from AMDA Re to SSCWeb R_E
  #info['xyz'] = xyz # Gives better match to CDAWeb
  info['time'] = hapitime2datetime(data['Time'])

  # Return not needed b/c info is modified in place. Keep for clarity.
  return info

def print_first_last_(sscweb_, cdaweb_, amda_, jpl_):
  if not print_first_last:
    return

  frame = sscweb_['frame']
  for i in [0, -1]:

    xyz_sscweb = sscweb_['xyz'][i]
    if cdaweb_ is not None:
      xyz_cdaweb = cdaweb_['xyz'][i]
    if jpl_ is not None:
      xyz_jpl = jpl_['xyz'][i]
    if amda_ is not None:
      xyz_amda = amda_['xyz'][i]

    time_sscweb = sscweb_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if cdaweb_ is not None:
      time_cdaweb = cdaweb_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if jpl_ is not None:
      time_jpl = jpl_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if amda_ is not None:
      time_amda = amda_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')

    u = 'R_E'
    if i == 0:
      print(f'  First values in {u}:')
    else:
      print(f'  Last values in {u}:')

    print('            {0:7s} {1:7s} {2:7s}'.format(f' X_{frame}', f' Y_{frame}', f' Z_{frame}'))
    print('    SSCWeb: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_sscweb, time_sscweb))
    if cdaweb_ is not None:
      print('    CDAWeb: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_cdaweb, time_cdaweb))
    if jpl_ is not None:
      print('       JPL: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_jpl, time_jpl))
    if amda_ is not None:
      print('      AMDA: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_amda, time_amda))

def interp(timei, time, xyz):
  from scipy.interpolate import interp1d

  # Interpolate CDAWeb data onto the union of timestamps
  interp_func = interp1d(
    [t.timestamp() for t in time],
    xyz,
    axis=0,
    bounds_error=False,
    fill_value="extrapolate"
  )

  return interp_func([t.timestamp() for t in timei])

def compute_diffs(info1, info2):

  print(f"  Computing differences between {info1['name']} and {info2['name']}")

  info1['time'] = [t.replace(tzinfo=None) for t in info1['time']]
  info2['time'] = [t.replace(tzinfo=None) for t in info2['time']]

  if interp_times:
    print("    Computing union of timestamps.")
    # Union of timestamps from CDAWeb and SSCWeb
    union_times = set(info1['time']).union(set(info2['time']))
    union_times = sorted(union_times)

    print(f"    Interpolating to {len(union_times)} timestamps.")
    info1['xyz'] = interp(union_times, info1['time'], info1['xyz'])
    info1['time'] = union_times
    info2['xyz'] = interp(union_times, info2['time'], info2['xyz'])
    info2['time'] = union_times
  else:
    # Common timestamps
    print("    Computing common timestamps.")
    common_times = set(info2['time']).intersection(set(info1['time']))
    common_times = sorted(common_times)
    common_times = sorted(common_times)
    print(f"    {len(common_times)} common timestamps.")

    # Find the indices of common times
    idx1 = np.nonzero(np.isin(info1['time'], common_times))[0]
    idx2 = np.nonzero(np.isin(info2['time'], common_times))[0]

    info1['xyz'] = np.array(info1['xyz'][idx1])
    info1['time'] = common_times
    info2['xyz'] = np.array(info2['xyz'][idx2])
    info2['time'] = common_times

  nt = len(info1['time'])
  Δθ = np.full(nt, np.nan)
  Δr = np.full(nt, np.nan)
  Δr_rel = np.full(nt, np.nan)
  t = np.empty(nt, dtype='datetime64[ns]')
  t[:] = np.datetime64('NaT')
  r = np.linalg.norm(info2['xyz'], axis=1)

  for i in range(nt):
    # Remove timezone from info_cdaweb['time']
    t[i] = np.datetime64(info1['time'][i])
    Δr[i] = np.linalg.norm(info1['xyz'][i] - info2['xyz'][i])

    Δr_rel[i] = Δr[i]/r[i]

    # d = denominator for angle calculation
    # Δθ[i] = arccos [ (a·b)/(|a|*|b|) ] = arccos (n/d)
    n = np.dot(info1['xyz'][i], info2['xyz'][i])
    d = np.linalg.norm(info1['xyz'][i])*np.linalg.norm(info2['xyz'][i])
    if np.abs(n) > np.abs(d):
      if np.abs(n) - np.abs(d) > 1e-10:
        raise ValueError(f"|n| - |d| > 1e-10. n = {n:.16f} d = {d:.16f}, n-d = {n-d:.16f}")
      if cos_warnings:
        # If n == d, then the angle is 0 degrees, but numpy gives
        # RuntimeWarning: invalid value encountered in arccos. Did rounding
        # cause n > d?
        # TODO: Consider https://github.com/sunpy/sunpy/pull/7530#issuecomment-2020890282
        wmsg =  f"    Warning: {t[i]}:\n"
        wmsg +=  "      |n| > |d| in Δθ[i] = arccos [ (a·b)/(|a|*|b|) ] = arccos (n/d)\n"
        wmsg += f"      a = {info1['xyz'][i]}\n"
        wmsg += f"      b = {info2['xyz'][i]}\n"
        wmsg += f"      n = {n:.16f}\n"
        wmsg += f"      d = {d:.16f}\n"
        wmsg += f"      n/d = {n/d}\n"
        wmsg +=  "      Using Δθ = 0\n"
        print(wmsg)
      Δθ[i] = 0
    else:
      Δθ[i] = (180/np.pi)*np.arccos(n/d)

    if print_interp_vals:
      t_str = info_cdaweb['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
      print(f"    t = {t_str} | Δr = {Δr[i]:.5f} [R_E] | Δ∠: {Δθ[i]:.5f}°")

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

def plot_xyz(ax, info1, info2):

  print(f"  Plotting xyz for {info1['name']} and {info2['name']}")

  colors = ['r', 'g', 'b']
  component = ['X', 'Y', 'Z']

  for c in range(3):
    label1 = f'{info1['name']}/{component[c]}'
    label2 = f'{info2['name']}/{component[c]}'
    ax.plot(info1['time'], info1['xyz'][:,c],
            label=label1, lw=2, linestyle='-', color=colors[c])
    ax.plot(info2['time'], info2['xyz'][:,c],
            label=label2, lw=3, linestyle='--', color=colors[c])

  r = np.linalg.norm(info1['xyz'], axis=1)
  label = '$\\overline{r}$'
  ax.plot(info1['time'], r, label=label, lw=2, linestyle='-', color='k')

  adjust_y_range(ax, gap_fraction=1)
  ax.set_ylabel('$R_E$', rotation=0)
  ax.grid()

  ax.legend(**{**legend_kwargs, 'ncol': 4})
  ax.set_xticklabels([])

def plot_diffs(ax, t, Δr, Δr_rel, Δθ):

  print("  Plotting diffs")

  lw = 2 # line width

  Δr_rel_max = np.nanmax(Δr_rel)
  #A = f'{Δr_rel_max:.1e}'.split('e')[0]
  #e = str(int(np.floor(np.log10(Δr_rel_max))))
  #Δr_rel_max_str = f" ($Δr_{{\\rm{{max}}}}/\\overline{{r}}$ = {A}·$10^{{{e}}}$)"
  Δr_rel_max_str = f" (max= 1/{1/Δr_rel_max:.0f})"

  Δr_max = np.nanmax(Δr)
  Δr_max_str = f'($Δr_{{\\rm{{max}}}} = {Δr_max*R_E:.1f}$ [km])'

  print(f"    Δr_max = {Δr_max:.5f} [R_E]")
  print(f"    Δr_max = {Δr_max*R_E:.1f} [km]")

  ax.plot(t, Δr, 'r-', lw=lw, label=f'$|Δr|/R_E$ {Δr_max_str}')
  ax.plot(t, Δr_rel, 'b-', lw=lw, label=f'$|Δr|/\\overline{{r}}$ {Δr_rel_max_str}')
  ax.plot(t, Δθ, 'g-', lw=lw, label='$Δθ$ [deg]')

  adjust_y_range(ax, bottom=0, gap_fraction=1)
  ax.legend(**legend_kwargs)
  ax.grid()

def _figprep():
  gs = plt.gcf().add_gridspec(2)
  axes = gs.subplots(sharex=True)

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

  return axes

def _savefigs(fname):
  for fmt in ['svg', 'png', 'pdf']:
    kwargs = {'bbox_inches': 'tight'}
    if fmt == 'png':
      kwargs['dpi'] = 300

    fname_full = f'figures/{fmt}/{fname}.{fmt}'
    os.makedirs(os.path.dirname(fname_full), exist_ok=True)
    print(f"  Writing {fname_full}")
    plt.savefig(fname_full, bbox_inches='tight')
  plt.close()

def _adjust_xlim(axes):
  from datetime import timedelta
  import matplotlib.dates as mdates
  from datetick import datetick

  datetick('x')#, debug=True)
  return
  # Subtract 1 minute to the x-axis limits to make sure the first tick is shown
  m = timedelta(minutes=1)
  xlim = [mdates.num2date(xlim) for xlim in axes[0].get_xlim()]
  hm = timedelta(hours=1, minutes=1)
  # Round to the next hour + 1 minute
  last = xlim[1].replace(minute=0, second=0, microsecond=0) + hm
  axes[0].set_xlim(xlim[0] - m, last)
  axes[1].set_xlim(xlim[0] - m, last)

  labels = axes[1].get_xticklabels()
  if labels:
    if "\n" in labels[0].get_text():
      labels[0].set_horizontalalignment('left')
    if "\n" in labels[-1].get_text():
      labels[-1].set_horizontalalignment('right')

def _plot(axes, satellite, frame, info1, info2):
  fname = f"{satellite}_{frame}_{info1['name']}_vs_{info2['name']}"
  axes[0].set_title(f"{satellite} {frame}")
  plot_xyz(axes[0], info1, info2)
  plot_diffs(axes[1], *compute_diffs(info1, info2))
  _adjust_xlim(axes)
  _savefigs(fname)

satellites = infos()

for satellite in satellites:

  if len(satellites_only) > 0 and satellite not in satellites_only:
    continue

  infos_ = infos(satellite)

  for i in range(len(infos_['sscweb'])):

    frame = infos_['sscweb'][i]['frame']

    print(f"\n{satellite} {frame}")

    # Adds xyz and time to info dict
    sscweb_ = sscweb(infos_['sscweb'][i], logging=hapi_logging)
    cdaweb_ = cdaweb(infos_['cdaweb'][i], logging=hapi_logging)
    amda_ = amda(infos_['amda'][i], logging=hapi_logging)
    jpl_ = jpl(infos_['sscweb'][i], logging=hapi_logging)

    print_first_last_(sscweb_, cdaweb_, amda_, jpl_)

    sources = [sscweb_, cdaweb_, amda_, jpl_]

    for i, source1 in enumerate(sources):
      for j, source2 in enumerate(sources):
        if i >= j or source1 is None or source2 is None:
          continue
        print(f"Comparing {source1['name']} with {source2['name']}")
        axes = _figprep()
        _plot(axes, satellite, frame, source1, source2)
