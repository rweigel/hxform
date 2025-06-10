# pip install netCDF4 'sunpy[all]' hapiclient

import datetime
import warnings
import sunpy
from netCDF4 import Dataset
from sunpy.coordinates import get_horizons_coord
from hapiclient import hapi
from hapiclient import hapitime2datetime

warnings.filterwarnings("ignore", message="The argument 'infer_datetime_format'")
hapi_logging = False
sunpy.log.setLevel('ERROR')

# NGDC data
# https://www.ngdc.noaa.gov/dscovr/portal/index.html#/download/1637802000000;1638666000000/pop
file_id = Dataset('data/nc/oe_pop_dscovr_s20211125000000_e20211125235959_p20211126022149_pub.nc', 'r')
print(file_id)
print(40 * '-')

for frame in ['GSE', 'J2K', 'GSM']:

  framelc = frame.lower()

  frame_cdaweb = frame
  frame_jpl = frame
  frame_ngdc = framelc
  frame_sscweb = frame
  if frame == 'J2K':
    frame_ngdc = 'gci'
    frame_cdaweb = 'GCI'

  time = file_id.variables['time'][:]
  x = file_id.variables[f'sat_x_{frame_ngdc}'][:]
  y = file_id.variables[f'sat_y_{frame_ngdc}'][:]
  z = file_id.variables[f'sat_z_{frame_ngdc}'][:]

  time_ngdc = datetime.datetime.fromtimestamp(time[0]/1000., tz=datetime.timezone.utc)
  time_ngdc = time_ngdc.strftime('%Y-%m-%dT%H:%M:%S.%fZ')
  xyz_ngdc = [x[0], y[0], z[0]]

  # JPL/Horizons web service
  if frame == 'GSE':
    q = get_horizons_coord('DSCOVR', time_ngdc).geocentricsolarecliptic.cartesian.xyz.to('km')
  if frame == 'J2K':
    q = get_horizons_coord('DSCOVR', time_ngdc).geocentricearthequatorial.cartesian.xyz.to('km')
  if frame == 'GSM':
    q = get_horizons_coord('DSCOVR', time_ngdc).geocentricsolarmagnetospheric.cartesian.xyz.to('km')
  xyz_jpl = q.value

  #solar_system_ephemeris.set('de432s')
  #>>> get_horizons_coord('DSCOVR', '2021-11-25').geocentricearthequatorial.cartesian.xyz.to('km')
  #<Quantity [ -705379.60527431, -1301735.49325424,  -473634.49321795] km>
  #>>> solar_system_ephemeris.set('de440s')
  #>>> get_horizons_coord('DSCOVR', '2021-11-25').geocentricearthequatorial.cartesian.xyz.to('km')
  # -705379.80      -1301735.36    -473634.52

  # SSCWeb via HAPI data
  # https://sscweb.gsfc.nasa.gov/cgi-bin/Locator.cgi?SPCR=dscovr&dscovr_res=&START_TIME=2021%2F11%2F25+0%3A0%3A0&STOP_TIME=2021%2F11%2F26+23%3A59%3A59&RESOLUTION=1&TOD=7&J2000=7&GEO=&GM=&GSE=&GSM=&SM=&REG_OPT=&MNMX_FLTR_ACCURACY=2&OPT=&TRC_GEON=&TRC_GEOS=&TRC_GMN=&TRC_GMS=&FILTER_DIST_UNITS=1&TOD_APPLY_FILTER=&TODX_MNMX=&TOD_XGT=&TOD_XLT=&TODY_MNMX=&TOD_YGT=&TOD_YLT=&TODZ_MNMX=&TOD_ZGT=&TOD_ZLT=&TODLAT_MNMX=&TOD_LATGT=&TOD_LATLT=&TODLON_MNMX=&TOD_LONGT=&TOD_LONLT=&TODLT_MNMX=&TOD_LTGT=&TOD_LTLT=&J2000_APPLY_FILTER=&J2000X_MNMX=&J2000_XGT=&J2000_XLT=&J2000Y_MNMX=&J2000_YGT=&J2000_YLT=&J2000Z_MNMX=&J2000_ZGT=&J2000_ZLT=&J2000LAT_MNMX=&J2000_LATGT=&J2000_LATLT=&J2000LON_MNMX=&J2000_LONGT=&J2000_LONLT=&J2000LT_MNMX=&J2000_LTGT=&J2000_LTLT=&GEO_APPLY_FILTER=&GEOX_MNMX=&GEO_XGT=&GEO_XLT=&GEOY_MNMX=&GEO_YGT=&GEO_YLT=&GEOZ_MNMX=&GEO_ZGT=&GEO_ZLT=&GEOLAT_MNMX=&GEO_LATGT=&GEO_LATLT=&GEOLON_MNMX=&GEO_LONGT=&GEO_LONLT=&GEOLT_MNMX=&GEO_LTGT=&GEO_LTLT=&GM_APPLY_FILTER=&GMX_MNMX=&GM_XGT=&GM_XLT=&GMY_MNMX=&GM_YGT=&GM_YLT=&GMZ_MNMX=&GM_ZGT=&GM_ZLT=&GMLAT_MNMX=&GM_LATGT=&GM_LATLT=&GMLON_MNMX=&GM_LONGT=&GM_LONLT=&GMLT_MNMX=&GM_LTGT=&GM_LTLT=&GSE_APPLY_FILTER=&GSEX_MNMX=&GSE_XGT=&GSE_XLT=&GSEY_MNMX=&GSE_YGT=&GSE_YLT=&GSEZ_MNMX=&GSE_ZGT=&GSE_ZLT=&GSELAT_MNMX=&GSE_LATGT=&GSE_LATLT=&GSELON_MNMX=&GSE_LONGT=&GSE_LONLT=&GSELT_MNMX=&GSE_LTGT=&GSE_LTLT=&GSM_APPLY_FILTER=&GSMX_MNMX=&GSM_XGT=&GSM_XLT=&GSMY_MNMX=&GSM_YGT=&GSM_YLT=&GSMZ_MNMX=&GSM_ZGT=&GSM_ZLT=&GSMLAT_MNMX=&GSM_LATGT=&GSM_LATLT=&GSMLON_MNMX=&GSM_LONGT=&GSM_LONLT=&GSMLT_MNMX=&GSM_LTGT=&GSM_LTLT=&SM_APPLY_FILTER=&SMX_MNMX=&SM_XGT=&SM_XLT=&SMY_MNMX=&SM_YGT=&SM_YLT=&SMZ_MNMX=&SM_ZGT=&SM_ZLT=&SMLAT_MNMX=&SM_LATGT=&SM_LATLT=&SMLON_MNMX=&SM_LONGT=&SM_LONLT=&SMLT_MNMX=&SM_LTGT=&SM_LTLT=&OTHER_FILTER_DIST_UNITS=1&RD_APPLY=&FS_APPLY=&NS_APPLY=&BS_APPLY=&MG_APPLY=&LV_APPLY=&IL_APPLY=&REG_FLTR_SWITCH=&SCR_APPLY=&SCR=&RTR_APPLY=&RTR=&BTR_APPLY=&NBTR=&SBTR=&EXTERNAL=3&EXT_T1989c=1&KP_LONG_89=4&INTERNAL=1&ALTITUDE=100&DAY=1&TIME=2&DISTANCE=1&DIST_DEC=2&DEG=1&DEG_DEC=8&DEG_DIR=1&OUTPUT_CDF=1&LINES_PAGE=1&RNG_FLTR_METHOD=&PREV_SECTION=FO&SSC=LOCATOR_GENERAL&SUBMIT=Submit+query+and+wait+for+output&.cgifields=DISTANCE&.cgifields=OUTPUT_CDF&.cgifields=DAY&.cgifields=TIME&.cgifields=DEG&.cgifields=DEG_DIR
  server     = 'http://hapi-server.org/servers/SSCWeb/hapi'
  dataset    = 'dscovr'

  parameters = f'X_{frame_sscweb},Y_{frame_sscweb},Z_{frame_sscweb}'

  start      = '2021-11-25T00:00:00Z'
  stop       = '2021-12-05T00:00:00Z'
  opts       = {'logging': hapi_logging, 'usecache': True, 'cachedir': './hapicache'}

  data_sscweb, meta_sscweb = hapi(server, dataset, parameters, start, stop, **opts)

  # Earth radius in km used by SSCWeb. See
  # https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
  # SSCWeb/HAPI provides ephemeris as three separate parameters. Here we combine parameters into a list.
  # After inspection of the metadata, we see that SSCWeb reports in R_E while CDAWeb in km.
  # SSCWeb reports in R_E
  # https://sscweb.gsfc.nasa.gov/sscweb_data_provenance.html
  R_E = 6378.16 # km
  xyz_sscweb = [R_E*data_sscweb[f'X_{frame_sscweb}'][1], R_E*data_sscweb[f'Y_{frame_sscweb}'][1], R_E*data_sscweb[f'Z_{frame_sscweb}'][1]]
  time_sscweb = hapitime2datetime([data_sscweb['Time'][0]])[0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')

  # CDAWeb via HAPI data
  if frame != 'GSM':
    server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
    dataset    = 'DSCOVR_ORBIT_PRE'
    parameters = f'{frame_cdaweb}_POS'
    start      = '2021-11-25T00:00:00Z'
    stop       = '2021-12-05T00:00:00Z'
    opts       = {'logging': hapi_logging, 'usecache': True, 'cachedir': './hapicache'}

    data_cdaweb, meta_cdaweb = hapi(server, dataset, parameters, start, stop, **opts)
    xyz_cdaweb = list(data_cdaweb[0][1])
    time_cdaweb = data_cdaweb['Time'][0].decode()

  print('        {0:13s} {1:13s} {2:13s}'.format(f'X_{frame} [km]', f'Y_{frame} [km]', f'Z_{frame} [km]'))
  print('SSCWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f}\t'.format(*xyz_sscweb, time_sscweb))
  if frame != 'GSM':
    print('CDAWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f}\t on {3:s}'.format(*xyz_cdaweb, time_cdaweb))
  print('NGDC:   {0:<13.3f} {1:<13.3f} {2:<13.3f}\t on {3:s}'.format(*xyz_ngdc, time_ngdc))
  print('JPL:    {0:<13.3f} {1:<13.3f} {2:<13.3f}\t on {3:s}'.format(*xyz_jpl, time_ngdc))
