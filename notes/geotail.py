#frame = 'GCI' # GCI, GSE
frame = 'GCI'

from hapiclient import hapi
#from hapiplot import hapiplot

server     = 'http://hapi-server.org/servers/SSCWeb/hapi'
dataset_sscweb = 'geotail'

if frame == 'GCI':
  parameters = f'X_TOD,Y_TOD,Z_TOD'
else:
  parameters = f'X_{frame},Y_{frame},Z_{frame}'

start      = '2021-11-25T00:00:00Z'
stop       = '2021-12-05T00:00:00Z'
opts       = {'logging': True, 'usecache': True, 'cachedir': './hapicache'}

data_sscweb1, meta_sscweb1 = hapi(server, dataset_sscweb, parameters, start, stop, **opts)
if frame == 'GCI':
  parameters = f'X_J2K,Y_J2K,Z_J2K'
  data_sscweb2, meta_sscweb2 = hapi(server, dataset_sscweb, parameters, start, stop, **opts)

# Earth radius in km used by SSCWeb; https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
R_E = 6378.16 
if frame == 'GCI':
  xyz_sscweb1 = [R_E*data_sscweb1[f'X_TOD'][1], R_E*data_sscweb1[f'Y_TOD'][1], R_E*data_sscweb1[f'Z_TOD'][1]]
  xyz_sscweb2 = [R_E*data_sscweb2[f'X_J2K'][1], R_E*data_sscweb2[f'Y_J2K'][1], R_E*data_sscweb2[f'Z_J2K'][1]]
else:
  xyz_sscweb1 = [R_E*data_sscweb1[f'X_{frame}'][1], R_E*data_sscweb1[f'Y_{frame}'][1], R_E*data_sscweb1[f'Z_{frame}'][1]]

server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
dataset    = 'GE_OR_DEF'
parameters = f'{frame}_POS'
start      = '2021-11-25T00:00:00Z'
stop       = '2021-12-05T00:00:00Z'
opts       = {'logging': True, 'usecache': True, 'cachedir': './hapicache'}

data_cdaweb, meta_cdaweb = hapi(server, dataset, parameters, start, stop, **opts)

from hapiclient import hapitime2datetime

xyz_cdaweb = list(data_cdaweb[0][1])

# THe HAPI SSCWeb server provides ephemeris as three separate parameters. Here we combine parameters into a list.
# After inspection of the metadata, we see that SSCWeb reports in R_E while CDAWeb in km.
# Convert CDAWeb data to R_E.

time_cdaweb = data_cdaweb['Time'][0].decode()

# Convert from YYYY-DOY to YYYY-MM-DD date format
time_sscweb = hapitime2datetime([data_sscweb1['Time'][0]])[0].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
# https://sscweb.gsfc.nasa.gov/sscweb_data_provenance.html
print('        {0:13s} {1:13s} {2:13s}'.format(f'X_{frame} [km]', f'Y_{frame} [km]', f'Z_{frame} [km]'))
print('CDAWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f} {3:s} (Using {4:s})'.format(*xyz_cdaweb, time_cdaweb, f'{dataset}/{frame}'))
if frame == 'GCI':
  print('SSCWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f} {3:s} (Using {4:s}/GEI/TOD)'.format(*xyz_sscweb1, time_sscweb, dataset_sscweb))
  print('SSCWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f} {3:s} (Using {4:s}/GEI/J2K)'.format(*xyz_sscweb2, time_sscweb, dataset_sscweb))
else:
  print('SSCWeb: {0:<13.3f} {1:<13.3f} {2:<13.3f} {3:s} (Using {4:s}/{5:s})'.format(*xyz_sscweb1, time_sscweb, dataset_sscweb, frame))

import numpy as np

xyz_cdaweb = np.array(xyz_cdaweb)
xyz_sscweb1 = np.array(xyz_sscweb1)
d = np.linalg.norm(xyz_cdaweb)*np.linalg.norm(xyz_sscweb1)
angle = (180/np.pi)*np.arccos(np.dot(xyz_cdaweb, xyz_sscweb1)/d)
print(angle)
#        X_GSE [R_E]  Y_GSE [R_E]  Z_GSE [R_E]
#SSCWeb: 243.36926843 1.43711310   13.03130715  	 on 2021-11-25T00:00:00.000000Z
#CDAWeb: 243.63903103 1.46400240   13.05533682  	 on 2021-11-25T00:00:00.000Z
