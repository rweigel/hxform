import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from hapiclient import hapi
from hapiclient import hapitime2datetime
import hapiplot

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# Earth radius in km used by SSCWeb. See
# https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
R_E = 6378.16 # km
satellite = 'mms'
satellite = 'geotail'

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

if satellite == 'mms':

  mms   = 'mms1'
  frame = 'GSE'
  #frame = 'GSM'

  start = '2016-09-01T00:00:00Z'
  stop  = '2016-09-02T00:00:00Z'

  info_cdaweb = {
    'dataset': 'MMS1_EPD-EIS_SRVY_L2_ELECTRONENERGY',
    'parameter': f'mms1_epd_eis_srvy_l2_electronenergy_position_{frame.lower()}',
    'start': start,
    'stop':  stop
  }

  info_sscweb = {
    'dataset': 'mms1',
    'frame': frame,
    'start': start,
    'stop':  stop
  }

if satellite == 'geotail':

  start      = '2021-11-25T00:00:00Z'
  stop       = '2021-12-05T00:00:00Z'

  infos_cdaweb = []
  infos_sscweb = []
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

  # https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html
  # https://cdaweb.gsfc.nasa.gov/misc/NotesG.html#GE_OR_DEF
  # CDAWeb dataset has
  #   GCI, GSE, GSM, HEC

def sscweb(info_sscweb):
  dataset    = info_sscweb['dataset']
  frame      = info_sscweb['frame']
  start      = info_sscweb['start']
  stop       = info_sscweb['stop']
  server = 'http://hapi-server.org/servers/SSCWeb/hapi'
  parameters = f'X_{frame},Y_{frame},Z_{frame}'

  opts       = {'logging': True, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameters, start, stop, **opts)

  xyz = np.column_stack((data[f'X_{frame}'], data[f'Y_{frame}'], data[f'Z_{frame}']))

  # Convert from YYYY-DOY to YYYY-MM-DD date format
  time = hapitime2datetime(data['Time'])

  return time, xyz

def cdaweb(info_cdaweb):
  dataset    = info_cdaweb['dataset']
  parameter  = info_cdaweb['parameter']
  start      = info_cdaweb['start']
  stop       = info_cdaweb['stop']
  server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
  opts       = {'logging': True, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameter, start, stop, **opts)

  xyz = data[parameter]
  xyz = xyz/R_E # Convert from km to R_E
  # THe HAPI SSCWeb server provides ephemeris as three separate parameters. Here
  # we combine parameters into a list. SSCWeb reports in R_E while CDAWeb in km.
  # Convert CDAWeb data to R_E.
  time = hapitime2datetime(data['Time'])

  return time, xyz

def print_first(time_cdaweb, time_sscweb, xyz_cdaweb, xyz_sscweb):
  time_cdaweb0 = time_cdaweb[0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  time_sscweb0 = time_sscweb[0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  # https://sscweb.gsfc.nasa.gov/sscweb_data_provenance.html
  u = '[R_E]'
  print('        {0:13s} {1:13s} {2:13s}'.format(f'X_{frame_cdaweb} {u}', f'Y_{frame_cdaweb} {u}', f'Z_{frame_cdaweb} {u}'))
  print('CDAWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f} {3:s} (Using {4:s}/{5:s})'.format(*xyz_cdaweb[0,:], time_cdaweb0, dataset_cdaweb, frame_cdaweb))
  print('SSCWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f} {3:s} (Using {4:s}/{5:s})'.format(*xyz_sscweb[0,:], time_sscweb0, dataset_sscweb, frame_sscweb))
  dr = np.linalg.norm(xyz_cdaweb[0,:] - xyz_sscweb[0,:])
  d = np.linalg.norm(xyz_cdaweb[0,:])*np.linalg.norm(xyz_sscweb[0,:])
  angle = (180/np.pi)*np.arccos(np.dot(xyz_cdaweb[0,:], xyz_sscweb[0,:])/d)
  #print(f"Δr = {dr/R_E:.3f} R_E = {dr:.0f} [km]; ∠: {angle:.3f}°")
  print(f"Δr = {dr:.3f} R_E; ∠: {angle:.3f}°")
  #        X_GSE [R_E]  Y_GSE [R_E]  Z_GSE [R_E]
  #SSCWeb: 243.36926843 1.43711310   13.03130715  	 on 2021-11-25T00:00:00.000000Z
  #CDAWeb: 243.63903103 1.46400240   13.05533682  	 on 2021-11-25T00:00:00.000Z

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

  interp_times = True
  if interp_times:
    # Create a union of timestamps from both CDAWeb and SSCWeb
    union_times = sorted(set(info_cdaweb['time']).union(set(info_sscweb['time'])))

    info_cdaweb['xyz'] = interp(union_times, info_cdaweb['time'], info_cdaweb['xyz'])
    info_cdaweb['time'] = union_times

    info_sscweb['xyz'] = interp(union_times, info_cdaweb['time'], info_sscweb['xyz'])
    info_sscweb['time'] = union_times

  # Find common timestamps
  common_times = set(info_cdaweb['time']).intersection(set(info_sscweb['time']))
  common_times = sorted(common_times)
  print(f"Common timestamps: {len(common_times)} found.")

  r_cdaweb = np.linalg.norm(info_cdaweb['xyz'], axis=1)
  r_diff = np.diff(r_cdaweb)/R_E
  t_diff = info_cdaweb['time'][1:]

  Δθ = np.full(len(common_times), np.nan)
  Δr = np.full(len(common_times), np.nan)
  Δr2 = np.full(len(common_times), np.nan)
  t = np.empty(len(common_times), dtype='datetime64[ns]')
  t[:] = np.datetime64('NaT')

  for i, tc in enumerate(common_times):

    # Find the index of common times
    idx_cdaweb = np.where(info_cdaweb['time'] == tc)[0][0]
    idx_sscweb = np.where(info_sscweb['time'] == tc)[0][0]

    xyz_cdaweb = np.array(info_cdaweb['xyz'][idx_cdaweb])
    xyz_sscweb = np.array(info_sscweb['xyz'][idx_sscweb])

    t[i] = np.datetime64(info_cdaweb['time'][idx_cdaweb])
    Δr[i] = np.linalg.norm(xyz_cdaweb - xyz_sscweb)#/R_E
    r2 = np.linalg.norm(xyz_cdaweb)
    Δr2[i] = np.linalg.norm(xyz_cdaweb - xyz_sscweb)/r2
    #Δr2[i] = np.linalg.norm(xyz_cdaweb)/R_E

    # d = denominator for angle calculation
    d = np.linalg.norm(xyz_cdaweb)*np.linalg.norm(xyz_sscweb)
    Δθ[i] = (180/np.pi)*np.arccos(np.dot(xyz_cdaweb, xyz_sscweb)/d)

    t_str = info_cdaweb['time'][idx_cdaweb].strftime('%Y-%m-%dT%H:%M:%S.%fZ')

    print(f"t = {t_str} | Δr = {Δr[i]:.5f} [R_E] | Δ∠: {Δθ[i]:.5f}°")

  return t, Δr, Δθ, t_diff, r_diff, Δr2

def plot_xyz(info_cdaweb, info_sscweb, title):
  colors = ['r', 'g', 'b']
  component = ['X', 'Y', 'Z']

  for c in range(3):
    plt.plot(info_cdaweb['time'], info_cdaweb['xyz'][:,c], label=f'CDAWeb/{component[c]}', lw=1, linestyle='-', color=colors[c])
    plt.plot(info_sscweb['time'], info_sscweb['xyz'][:,c], label=f'SSCWeb/{component[c]}', lw=2, linestyle='--', color=colors[c])

  plt.title(title, fontsize=12)
  plt.grid()
  plt.legend()
  hapiplot.plot.datetick.datetick(dir='x')

  fname = f'figures/{info_sscweb['dataset']}_{info_sscweb['frame']}'
  print(f"Writing {fname}.{{svg,png}}")
  plt.savefig(f'{fname}.svg', bbox_inches='tight')
  plt.savefig(f'{fname}.png', dpi=300, bbox_inches='tight')
  plt.close()

def plot_diffs(info_cdaweb, info_sscweb, title):
  t, Δr, Δθ, t_diff, r_diff, Δr2 = compute_diffs(info_cdaweb, info_sscweb)

  #bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1')
  bbox = None
  ms = 2
  plt.title(title, fontsize=12)

  plt.plot(t, Δr, 'k-', ms=ms, label='$Δr/R_E$')
  #plt.text(t[-1], Δr[-1], '  $Δr/R_E$', color='k', fontsize=12, bbox=bbox)

  plt.plot(t, Δθ, 'g-', ms=ms, label='$Δθ$ [deg]')
  #plt.text(t[-1], Δθ[-1], '  $Δθ$ [deg]', color='g', fontsize=12, bbox=bbox)

  #plt.plot(t_diff, np.abs(r_diff), 'r.', ms=ms, label='$|r_{t+1}-r_t|/R_E$')
  #plt.text(t_diff[-1], np.abs(r_diff[-1]), '  $|r_{t+1}-r_t|/R_E$', color='r', fontsize=12, bbox=bbox, ha='left', va='top')

  plt.plot(t, 10*Δr2, 'b-', ms=ms, label='$10Δr/r$')
  #plt.text(t[-1], 10*Δr2[-1], '  $10Δr/r$', color='b', fontsize=12, bbox=bbox)

  #ax = plt.gca()
  #ax.spines['top'].set_visible(False)
  #ax.spines['right'].set_visible(False)
  plt.legend(loc='upper left', markerscale=3, fontsize=12, scatterpoints=3, frameon=True, ncol=3)
  plt.grid()
  hapiplot.plot.datetick.datetick(dir='x')

  ylim = plt.gca().get_ylim()
  yticks = plt.gca().get_yticks()
  half_width = (yticks[1] - yticks[0]) / 2
  plt.ylim(0, ylim[1] + half_width)

  fname = f'figures/{info_sscweb['dataset']}_{info_sscweb['frame']}-diffs'
  print(f"Writing {fname}.{{svg,png}}")
  plt.savefig(f'{fname}.svg', bbox_inches='tight')
  plt.savefig(f'{fname}.png', dpi=300, bbox_inches='tight')


for i in range(len(infos_cdaweb)):
  info_cdaweb = infos_cdaweb[i]
  info_sscweb = infos_sscweb[i]

  info_sscweb['time'], info_sscweb['xyz'] = sscweb(info_sscweb)
  info_cdaweb['time'], info_cdaweb['xyz'] = cdaweb(info_cdaweb)

  title = f"SSCWeb/{info_sscweb['dataset']}/{info_sscweb['frame']} | "
  title += f"CDAWeb/{info_cdaweb['dataset']}/{info_cdaweb['parameter']}"

  plot_xyz(info_cdaweb, info_sscweb, title)
  plot_diffs(info_cdaweb, info_sscweb, title)
