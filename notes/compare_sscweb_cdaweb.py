import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from hapiclient import hapi
from hapiclient import hapitime2datetime
import hapiplot

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.figsize'] = (8.5, 11)

legend_kwargs = {
  'loc': 'upper center',
  'markerscale': 3,
  'borderaxespad': 0,
  'framealpha': 1,
  'fontsize': 12,
  'frameon': True,
  'ncol': 3
}

# Earth radius in km used by SSCWeb. See
# https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
R_E = 6378.16 # km
satellite = '' # Set to '' to run all
#satellite = 'mms'
#satellite = 'geotail'
satellite = 'themis'
hapi_logging = True
interp_times = True
print_vals = False
print_first = True

"""
SSCWeb provides TOD, J2K, GEO, GM, GSE, and SM
  TOD in SSCWeb HAPI server is same as GEI from SSCWeb.
  (The reason the HAPI server uses TOD for GEI is that "TOD" is the 
  POST query parameter names used to request data in GEI is "TOD" and not GEI
  and consistency in the query parameters was sought).

  From https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html, 
  "Geocentric Equatorial Inertial system. This system has X-axis
  pointing from the Earth toward the first point of Aries (the position of
  the Sun at the vernal equinox). This direction is the intersection of the
  Earth's equatorial plane and the ecliptic plane and thus the X-axis lies
  in both planes. The Z-axis is parallel to the rotation axis of the Earth,
  and y completes the right-handed orthogonal set (Y = Z * X). Geocentric
  Inertial (GCI) and Earth-Centered Inertial (ECI) are the same as GEI."

  Based on the above quote, we treat GCI from CDAWeb as the same as GEI from SSCWeb.
"""

infos_cdaweb = []
infos_sscweb = []

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
      frame = 'TOD'

    infos_sscweb.append({
      'dataset': themis,
      'frame': frame,
      'start': start,
      'stop':  stop
    })

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

  # https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html
  # https://cdaweb.gsfc.nasa.gov/misc/NotesG.html#GE_OR_DEF
  # CDAWeb dataset has
  #   GCI, GSE, GSM, HEC

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
      frame = 'TOD'

    infos_sscweb.append({
      'dataset': 'geotail',
      'frame': frame,
      'start': start,
      'stop':  stop
    })

def sscweb(info_sscweb, logging=False):
  dataset    = info_sscweb['dataset']
  frame      = info_sscweb['frame']
  start      = info_sscweb['start']
  stop       = info_sscweb['stop']
  server = 'http://hapi-server.org/servers/SSCWeb/hapi'
  parameters = f'X_{frame},Y_{frame},Z_{frame}'

  opts       = {'logging': logging, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameters, start, stop, **opts)

  xyz = np.column_stack((data[f'X_{frame}'], data[f'Y_{frame}'], data[f'Z_{frame}']))

  # Convert from YYYY-DOY to YYYY-MM-DD date format
  time = hapitime2datetime(data['Time'])

  return time, xyz

def cdaweb(info_cdaweb, logging=False):
  dataset    = info_cdaweb['dataset']
  parameter  = info_cdaweb['parameter']
  start      = info_cdaweb['start']
  stop       = info_cdaweb['stop']
  server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
  opts       = {'logging': logging, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameter, start, stop, **opts)

  xyz = data[parameter]
  xyz = xyz/R_E # Convert from km to R_E
  # THe HAPI SSCWeb server provides ephemeris as three separate parameters. Here
  # we combine parameters into a list. SSCWeb reports in R_E while CDAWeb in km.
  # Convert CDAWeb data to R_E.
  time = hapitime2datetime(data['Time'])

  return time, xyz

def print_first_(info_cdaweb, info_sscweb):
  if not print_first:
    return

  frame = info_sscweb['frame']
  xyz_cdaweb0 = info_cdaweb['xyz'][0]
  xyz_sscweb0 = info_sscweb['xyz'][0]
  time_cdaweb0 = info_cdaweb['time'][0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  time_sscweb0 = info_sscweb['time'][0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')

  u = 'R_E'
  print(f'  First values in {u}:')
  print('            {0:7s} {1:7s} {2:7s}'.format(f'X_{frame}', f'Y_{frame}', f'Z_{frame}'))
  print('    CDAWeb: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_cdaweb0, time_cdaweb0))
  print('    SSCWeb: {0:<7.3f} {1:<7.3f} {2:<7.3f} {3:s}'.format(*xyz_sscweb0, time_sscweb0))

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
  return (ra + rb)/2

def plot_xyz(ax, info_cdaweb, info_sscweb):
  colors = ['r', 'g', 'b']
  component = ['X', 'Y', 'Z']

  for c in range(3):
    ax.plot(info_cdaweb['time'], info_cdaweb['xyz'][:,c], label=f'CDAWeb/{component[c]}', lw=2, linestyle='-', color=colors[c])
    ax.plot(info_sscweb['time'], info_sscweb['xyz'][:,c], label=f'SSCWeb/{component[c]}', lw=3, linestyle='--', color=colors[c])

  r = compute_r(info_cdaweb, info_sscweb)
  ax.plot(info_cdaweb['time'], r, label='$\\overline{r}$', lw=2, linestyle='-', color='k')

  adjust_y_range(ax, gap_fraction=0.75)
  ax.set_ylabel('$R_E$', fontsize=12, rotation=0)
  ax.grid()
  ax.legend(**{**legend_kwargs, 'ncol': 4})
  #ax.legend(**{**legend_kwargs, 'ncol': 4})
  ax.set_xticklabels([])

def plot_diffs(ax, info_cdaweb, info_sscweb):
  print("  Computing differences")
  t, Δr, Δr_rel, Δθ  = compute_diffs(info_cdaweb, info_sscweb)

  lw = 2 # line width

  Δr_rel_max = np.nanmax(Δr_rel)
  print(f"  Δr_max = {Δr_rel_max:.5f} [R_E]")
  A = f'{Δr_rel_max:.1e}'.split('e')[0]
  Δr_rel_max_str = f"{A}·10^{{{int(np.floor(np.log10(Δr_rel_max)))}}}"  # Convert to 10^ notation
  ax.plot(t, Δr, 'r-', lw=lw, label='$Δr/R_E$')
  ax.plot(t, Δθ, 'g-', lw=lw, label='$Δθ$ [deg]')
  ax.plot(t, Δr_rel, 'b-', lw=lw, label=f'$Δr/\\overline{{r}}$  (max = ${Δr_rel_max_str}$)')

  adjust_y_range(ax, bottom=0)
  ax.legend(**legend_kwargs)
  ax.grid()

for i in range(len(infos_cdaweb)):

  info_cdaweb = infos_cdaweb[i]
  info_sscweb = infos_sscweb[i]

  info_cdaweb['dataset']
  a = f"SSCWeb/{info_sscweb['dataset']}/{info_sscweb['frame']}"
  b = f"CDAWeb/{info_cdaweb['dataset'].split('@')[0]}/{info_cdaweb['parameter']}"
  if info_sscweb['dataset'].startswith('mms'):
    # Shorten name for MMS
    b = f"CDAWeb/{info_cdaweb['parameter']}"
  print(f"Comparing\n{a}\nwith\n{b}")
  title = f"{a} vs. {b}"

  info_sscweb['time'], info_sscweb['xyz'] = sscweb(info_sscweb, logging=hapi_logging)
  info_cdaweb['time'], info_cdaweb['xyz'] = cdaweb(info_cdaweb, logging=hapi_logging)

  print_first_(info_cdaweb, info_sscweb)

  gs = plt.gcf().add_gridspec(2)
  axes = gs.subplots(sharex=True)
  axes[0].set_title(title, fontsize=12)

  for ax in axes:
    ax.spines['bottom'].set_visible(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', which='minor', length=0)  # Remove minor tick lines next to y-axis numbers
    ax.tick_params(axis='y', length=0)  # Remove tick lines next to y-axis numbers

    for ax in axes:
      ax.grid(which='minor', linestyle=':', linewidth=0.5, color=[0.75]*3)
      ax.minorticks_on()
  print("  Plotting xyz")
  plot_xyz(axes[0], info_cdaweb, info_sscweb)
  print("  Plotting difs")
  plot_diffs(axes[1], info_cdaweb, info_sscweb)

  from datetime import timedelta
  # Add 1 minute to the x-axis limits to make sure the last tick is shown
  axes[0].set_xlim(info_cdaweb['time'][0], info_cdaweb['time'][-1] + timedelta(minutes=1))
  axes[1].set_xlim(info_cdaweb['time'][0], info_cdaweb['time'][-1] + timedelta(minutes=1))
  hapiplot.plot.datetick.datetick(dir='x')

  labels = axes[1].get_xticklabels()
  if labels:
    labels[0].set_horizontalalignment('left')
    labels[-1].set_horizontalalignment('right')

  fname = f'figures/{info_sscweb['dataset']}_{info_sscweb['frame']}'
  #plt.savefig(f'{fname}.svg', bbox_inches='tight')
  #plt.savefig(f'{fname}.png', dpi=300, bbox_inches='tight')
  print(f"  Writing {fname}.pdf")
  plt.savefig(f'{fname}.pdf', bbox_inches='tight')
  plt.close()

  print("")
