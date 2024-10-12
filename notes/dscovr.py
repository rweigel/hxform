#frame = 'GCI' # GCI, GSE
frame = 'GSE'
framelc = frame.lower()

# Library to work with netCDF files
from netCDF4 import Dataset

# Open a .nc file ("file_name")
# https://www.ngdc.noaa.gov/dscovr/portal/index.html#/download/1637802000000;1638666000000/pop
file_id = Dataset('nc/oe_pop_dscovr_s20211125000000_e20211125235959_p20211126022149_pub.nc', 'r')
print(file_id)
time = file_id.variables['time'][:]
x = file_id.variables[f'sat_x_{framelc}'][:]
y = file_id.variables[f'sat_y_{framelc}'][:]
z = file_id.variables[f'sat_z_{framelc}'][:]

from datetime import datetime

time_ngdc = datetime.utcfromtimestamp(time[0]/1000.).strftime('%Y-%m-%dT%H:%M:%S.%fZ')
xyz_ngdc = [x[0], y[0], z[0]]

#print(time_ngdc, xyz_ngdc)
#exit()
#for i in range(len(time)):
#  print(datetime.utcfromtimestamp(time[i]/1000.), x[i], y[i], z[i])

#from sunpy.coordinates import get_horizons_coord
#>>> get_horizons_coord('DSCOVR', '2021-11-25').geocentricearthequatorial.cartesian.xyz.to('km')
#<Quantity [ -705386.15221142, -1301734.17689259,  -473629.32913673] km>

#from astropy.coordinates import solar_system_ephemeris
#>>> solar_system_ephemeris.set('de432s')
#>>> get_horizons_coord('DSCOVR', '2021-11-25').geocentricearthequatorial.cartesian.xyz.to('km')
#<Quantity [ -705379.60527431, -1301735.49325424,  -473634.49321795] km>

# DE440s
#>>> solar_system_ephemeris.set('de440s')
#>>> get_horizons_coord('DSCOVR', '2021-11-25').geocentricearthequatorial.cartesian.xyz.to('km')
# -705379.80      -1301735.36    -473634.52

# 01a
from hapiclient import hapi
#from hapiplot import hapiplot

# https://sscweb.gsfc.nasa.gov/cgi-bin/Locator.cgi?SPCR=dscovr&dscovr_res=&START_TIME=2021%2F11%2F25+0%3A0%3A0&STOP_TIME=2021%2F11%2F26+23%3A59%3A59&RESOLUTION=1&TOD=7&J2000=7&GEO=&GM=&GSE=&GSM=&SM=&REG_OPT=&MNMX_FLTR_ACCURACY=2&OPT=&TRC_GEON=&TRC_GEOS=&TRC_GMN=&TRC_GMS=&FILTER_DIST_UNITS=1&TOD_APPLY_FILTER=&TODX_MNMX=&TOD_XGT=&TOD_XLT=&TODY_MNMX=&TOD_YGT=&TOD_YLT=&TODZ_MNMX=&TOD_ZGT=&TOD_ZLT=&TODLAT_MNMX=&TOD_LATGT=&TOD_LATLT=&TODLON_MNMX=&TOD_LONGT=&TOD_LONLT=&TODLT_MNMX=&TOD_LTGT=&TOD_LTLT=&J2000_APPLY_FILTER=&J2000X_MNMX=&J2000_XGT=&J2000_XLT=&J2000Y_MNMX=&J2000_YGT=&J2000_YLT=&J2000Z_MNMX=&J2000_ZGT=&J2000_ZLT=&J2000LAT_MNMX=&J2000_LATGT=&J2000_LATLT=&J2000LON_MNMX=&J2000_LONGT=&J2000_LONLT=&J2000LT_MNMX=&J2000_LTGT=&J2000_LTLT=&GEO_APPLY_FILTER=&GEOX_MNMX=&GEO_XGT=&GEO_XLT=&GEOY_MNMX=&GEO_YGT=&GEO_YLT=&GEOZ_MNMX=&GEO_ZGT=&GEO_ZLT=&GEOLAT_MNMX=&GEO_LATGT=&GEO_LATLT=&GEOLON_MNMX=&GEO_LONGT=&GEO_LONLT=&GEOLT_MNMX=&GEO_LTGT=&GEO_LTLT=&GM_APPLY_FILTER=&GMX_MNMX=&GM_XGT=&GM_XLT=&GMY_MNMX=&GM_YGT=&GM_YLT=&GMZ_MNMX=&GM_ZGT=&GM_ZLT=&GMLAT_MNMX=&GM_LATGT=&GM_LATLT=&GMLON_MNMX=&GM_LONGT=&GM_LONLT=&GMLT_MNMX=&GM_LTGT=&GM_LTLT=&GSE_APPLY_FILTER=&GSEX_MNMX=&GSE_XGT=&GSE_XLT=&GSEY_MNMX=&GSE_YGT=&GSE_YLT=&GSEZ_MNMX=&GSE_ZGT=&GSE_ZLT=&GSELAT_MNMX=&GSE_LATGT=&GSE_LATLT=&GSELON_MNMX=&GSE_LONGT=&GSE_LONLT=&GSELT_MNMX=&GSE_LTGT=&GSE_LTLT=&GSM_APPLY_FILTER=&GSMX_MNMX=&GSM_XGT=&GSM_XLT=&GSMY_MNMX=&GSM_YGT=&GSM_YLT=&GSMZ_MNMX=&GSM_ZGT=&GSM_ZLT=&GSMLAT_MNMX=&GSM_LATGT=&GSM_LATLT=&GSMLON_MNMX=&GSM_LONGT=&GSM_LONLT=&GSMLT_MNMX=&GSM_LTGT=&GSM_LTLT=&SM_APPLY_FILTER=&SMX_MNMX=&SM_XGT=&SM_XLT=&SMY_MNMX=&SM_YGT=&SM_YLT=&SMZ_MNMX=&SM_ZGT=&SM_ZLT=&SMLAT_MNMX=&SM_LATGT=&SM_LATLT=&SMLON_MNMX=&SM_LONGT=&SM_LONLT=&SMLT_MNMX=&SM_LTGT=&SM_LTLT=&OTHER_FILTER_DIST_UNITS=1&RD_APPLY=&FS_APPLY=&NS_APPLY=&BS_APPLY=&MG_APPLY=&LV_APPLY=&IL_APPLY=&REG_FLTR_SWITCH=&SCR_APPLY=&SCR=&RTR_APPLY=&RTR=&BTR_APPLY=&NBTR=&SBTR=&EXTERNAL=3&EXT_T1989c=1&KP_LONG_89=4&INTERNAL=1&ALTITUDE=100&DAY=1&TIME=2&DISTANCE=1&DIST_DEC=2&DEG=1&DEG_DEC=8&DEG_DIR=1&OUTPUT_CDF=1&LINES_PAGE=1&RNG_FLTR_METHOD=&PREV_SECTION=FO&SSC=LOCATOR_GENERAL&SUBMIT=Submit+query+and+wait+for+output&.cgifields=DISTANCE&.cgifields=OUTPUT_CDF&.cgifields=DAY&.cgifields=TIME&.cgifields=DEG&.cgifields=DEG_DIR
server     = 'http://hapi-server.org/servers/SSCWeb/hapi'
dataset    = 'dscovr'

if frame == 'GCI':
  parameters = f'X_TOD,Y_TOD,Z_TOD'
  #parameters = f'X_J2K,Y_J2K,Z_J2K'
else:
  parameters = f'X_{frame},Y_{frame},Z_{frame}'

start      = '2021-11-25T00:00:00Z'
stop       = '2021-12-05T00:00:00Z'
opts       = {'logging': True, 'usecache': True, 'cachedir': './hapicache'}

data_sscweb, meta_sscweb = hapi(server, dataset, parameters, start, stop, **opts)

R_E = 6378.16 # Earth radius in km
if frame == 'GCI':
  xyz_sscweb = [R_E*data_sscweb[f'X_TOD'][1], R_E*data_sscweb[f'Y_TOD'][1], R_E*data_sscweb[f'Z_TOD'][1]]
  #xyz_sscweb = [R_E*data_sscweb[f'X_J2K'][1], R_E*data_sscweb[f'Y_J2K'][1], R_E*data_sscweb[f'Z_J2K'][1]]
else:
  xyz_sscweb = [R_E*data_sscweb[f'X_{frame}'][1], R_E*data_sscweb[f'Y_{frame}'][1], R_E*data_sscweb[f'Z_{frame}'][1]]

server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
dataset    = 'DSCOVR_ORBIT_PRE'
parameters = f'{frame}_POS'
start      = '2021-11-25T00:00:00Z'
stop       = '2021-12-05T00:00:00Z'
opts       = {'logging': True, 'usecache': True, 'cachedir': './hapicache'}

data_cdaweb, meta_cdaweb = hapi(server, dataset, parameters, start, stop, **opts)

from hapiclient import hapitime2datetime

xyz_cdaweb = list(data_cdaweb[0][1])

# SSCWeb provides ephemeris as three separate parameters. Here we combine parameters into a list.
# After inspection of the metadata, we see that SSCWeb reports in R_E while CDAWeb in km.
# Convert CDAWeb data to R_E.

time_cdaweb = data_cdaweb['Time'][0].decode()

# Convert from YYYY-DOY to YYYY-MM-DD date format
time_sscweb = hapitime2datetime([data_sscweb['Time'][0]])[0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
# https://sscweb.gsfc.nasa.gov/sscweb_data_provenance.html
print('        {0:13s} {1:13s} {2:13s}'.format(f'X_{frame} [km]', f'Y_{frame} [km]', f'Z_{frame} [km]'))
print('SSCWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f}\t on {3:s}'.format(*xyz_sscweb, time_sscweb))
print('CDAWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f}\t on {3:s}'.format(*xyz_cdaweb, time_cdaweb))
print('NGDC:   {0:<13.3f} {1:<13.3f} {2:<13.3f}\t on {3:s}'.format(*xyz_ngdc, time_ngdc))

#        X_GSE [R_E]  Y_GSE [R_E]  Z_GSE [R_E]
#SSCWeb: 243.36926843 1.43711310   13.03130715  	 on 2021-11-25T00:00:00.000000Z
#CDAWeb: 243.63903103 1.46400240   13.05533682  	 on 2021-11-25T00:00:00.000Z


#hapiplot(data_cdaweb, meta_cdaweb);
#hapiplot(data_sscweb, meta_sscweb);



