import numpy as np
from hapiclient import hapi
from hapiclient import hapitime2datetime
import hapiplot

if True:
  frame_cdaweb = 'GSE'
  frame_sscweb = 'GSE'
  frame_cdaweb = 'GSM'
  frame_sscweb = 'GSM'

  #dataset_sscweb = 'mms3'
  #dataset_cdaweb = 'MMS3_EPD-EIS_SRVY_L2_ELECTRONENERGY'
  #parameter_cdaweb = 'mms3_epd_eis_srvy_l2_electronenergy_position_gse'
  #dataset_sscweb = 'mms2'
  #dataset_cdaweb = 'MMS2_EPD-EIS_SRVY_L2_ELECTRONENERGY'
  #parameter_cdaweb = 'mms2_epd_eis_srvy_l2_electronenergy_position_gse'
  dataset_sscweb = 'mms1'
  dataset_cdaweb = 'MMS1_EPD-EIS_SRVY_L2_ELECTRONENERGY'
  parameter_cdaweb = f'mms1_epd_eis_srvy_l2_electronenergy_position_{frame_cdaweb.lower()}'

  start      = '2016-09-01T00:00:00Z'
  stop       = '2016-09-02T00:00:00Z'

if False:
  dataset_sscweb = 'geotail'
  # https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html
  dataset_cdaweb = 'GE_OR_DEF'
  # https://cdaweb.gsfc.nasa.gov/misc/NotesG.html#GE_OR_DEF
  # CDAWeb dataset has
  #   GCI, GSE, GSM, HEC

  start      = '2021-11-25T00:00:00Z'
  stop       = '2021-12-05T00:00:00Z'

  #frame_cdaweb = 'GSM'
  #frame_sscweb = 'GSM'

  #frame_cdaweb = 'GSE'
  #frame_sscweb = 'GSE'

  #frame_cdaweb = 'GCI'
  #frame_sscweb = 'TOD'

  parameter_cdaweb = f'{frame_cdaweb}_POS'

# Earth radius in km used by SSCWeb. See
# https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
R_E = 6378.16 # km

# SSCWeb provides
#   TOD, J2K, GEO, GM, GSE, SM
#   TOD in SSCWeb HAPI server is same as GEI from SSCWeb.
#   (The reason the HAPI server uses TOD for GEI is that "TOD" is the
#   used the POST query parameters to request data in GEI and consistency in the
#   query parameters was sought).
#   From https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html, 
#   "Geocentric Equatorial Inertial system. This system has X-axis
#   pointing from the Earth toward the first point of Aries (the position of
#   the Sun at the vernal equinox). This direction is the intersection of the
#   Earth's equatorial plane and the ecliptic plane and thus the X-axis lies
#   in both planes. The Z-axis is parallel to the rotation axis of the Earth,
#   and y completes the right-handed orthogonal set (Y = Z * X). Geocentric
#   Inertial (GCI) and Earth-Centered Inertial (ECI) are the same as GEI."
#
#   Based on the above quote, we treat GCI from CDAWeb as the same as GEI from SSCWeb.

def sscweb(dataset, frame, start, stop):
  server = 'http://hapi-server.org/servers/SSCWeb/hapi'
  parameters = f'X_{frame},Y_{frame},Z_{frame}'

  opts       = {'logging': True, 'usecache': True, 'cachedir': './hapicache'}

  data, meta = hapi(server, dataset, parameters, start, stop, **opts)

  xyz = np.column_stack((data[f'X_{frame}'], data[f'Y_{frame}'], data[f'Z_{frame}']))

  # Convert from YYYY-DOY to YYYY-MM-DD date format
  time = hapitime2datetime(data['Time'])

  return time, xyz

def cdaweb(dataset, parameter, start, stop):
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

time_sscweb, xyz_sscweb = sscweb(dataset_sscweb, frame_sscweb, start, stop)
time_cdaweb, xyz_cdaweb = cdaweb(dataset_cdaweb, parameter_cdaweb, start, stop)

from matplotlib import pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
title = f'SSCWeb/{dataset_sscweb}/{frame_sscweb} | CDAWeb/{dataset_cdaweb}/{frame_cdaweb}'
colors = ['r', 'g', 'b']
component = ['X', 'Y', 'Z']
for c in range(3):
  plt.plot(time_cdaweb, xyz_cdaweb[:,c], label=f'CDAWeb/{component[c]}', lw=1, linestyle='-', color=colors[c])
  plt.plot(time_sscweb, xyz_sscweb[:,c], label=f'SSCWeb/{component[c]}', lw=2, linestyle='--', color=colors[c])
plt.title(title, fontsize=12)
plt.grid()
plt.legend()
hapiplot.plot.datetick.datetick(dir='x')
fname = f'figures/{dataset_sscweb}_{frame_sscweb}'
print(f"Writing {fname}.{{svg,png}}")
plt.savefig(f'{fname}.svg', bbox_inches='tight')
plt.savefig(f'{fname}.png', dpi=300, bbox_inches='tight')
plt.close()

if True:
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

exit()

# Find timestamps in time_cdaweb that are in time_sscweb
common_times = set(time_cdaweb).intersection(set(time_sscweb))
print(f"Common timestamps: {len(common_times)} found.")
common_times = sorted(common_times)

def compute_diffs(time_cdaweb, time_sscweb, xyz_cdaweb, xyz_sscweb):
  # Find timestamps in time_cdaweb that are in time_sscweb
  common_times = set(time_cdaweb).intersection(set(time_sscweb))
  print(f"Common timestamps: {len(common_times)} found.")
  common_times = sorted(common_times)

  r_cdaweb = np.linalg.norm(xyz_cdaweb, axis=1)
  r_diff = np.diff(r_cdaweb)/R_E
  t_diff = time_cdaweb[1:]

  Δθ = np.full(len(common_times), np.nan)
  Δr = np.full(len(common_times), np.nan)
  Δr2 = np.full(len(common_times), np.nan)
  t = np.empty(len(common_times), dtype='datetime64[ns]')
  t[:] = np.datetime64('NaT')

  for i, tc in enumerate(common_times):

    # Find the index of the timestamp in time_cdaweb
    idx_cdaweb = np.where(time_cdaweb == tc)[0][0]
    # Find the index of the timestamp in time_sscweb
    idx_sscweb = np.where(time_sscweb == tc)[0][0]

    xyz_cdaweb_t = np.array(xyz_cdaweb[idx_cdaweb])
    xyz_sscweb_t = np.array(xyz_sscweb[idx_sscweb])

    t[i] = np.datetime64(time_cdaweb[idx_cdaweb])
    Δr[i] = np.linalg.norm(xyz_cdaweb_t - xyz_sscweb_t)#/R_E
    r2 = np.linalg.norm(xyz_cdaweb_t)
    Δr2[i] = np.linalg.norm(xyz_cdaweb_t - xyz_sscweb_t)/r2
    #Δr2[i] = np.linalg.norm(xyz_cdaweb_t)/R_E

    # d = denominator for angle calculation
    d = np.linalg.norm(xyz_cdaweb_t)*np.linalg.norm(xyz_sscweb_t)
    Δθ[i] = (180/np.pi)*np.arccos(np.dot(xyz_cdaweb_t, xyz_sscweb_t)/d)

    t_str = time_cdaweb[idx_cdaweb].strftime('%Y-%m-%dT%H:%M:%S.%fZ')

    print(f"t = {t_str} Δr = {Δr[i]:.3f} [R_E] ∠: {Δθ[i]:.3f}°")

  return t, Δr, Δθ, t_diff, r_diff, Δr2

t, Δr, Δθ, t_diff, r_diff, Δr2 = compute_diffs(time_cdaweb, time_sscweb, xyz_cdaweb, xyz_sscweb)


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

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(loc='upper left', markerscale=3, fontsize=12, scatterpoints=3, frameon=True, ncol=3)
plt.grid()
hapiplot.plot.datetick.datetick(dir='x')
fname = f'figures/{dataset_sscweb}_{frame_sscweb}'
print(f"Writing {fname}.{{svg,png}}")
plt.savefig(f'{fname}-diffs.svg', bbox_inches='tight')
plt.savefig(f'{fname}-diffs.png', dpi=300, bbox_inches='tight')