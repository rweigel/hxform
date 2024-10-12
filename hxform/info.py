def lib_info(lib):
  infos = known_libs(info=True)
  for info in infos:
    if info['name'] == lib:
      return info
  return None


def known_transforms(lib):
  return lib_info(lib)["systems"]


def known_libs(info=False):

  import sunpy
  import spacepy
  import spiceypy

  # https://spacepy.github.io/coordinates.html
  systems_spacepy = ["ECI2000", "ECIMOD", "ECITOD", "GEI", "GSE", "GSM", "GEO", "SM", "MAG"] 

  knowns = [
            {
              'name': 'cxform',
              'version': None, # TODO: Use commit hash.
              'version_info': None,
              'systems': ["J2000", "GEI", "GEO", "GSE", "GSM", "SM", "MAG", "RTN", "GSEQ", "HEE", "HAE", "HEEQ"]
            },
            {
              'name': 'geopack_08_dp',
              'version': None,
              'version_info': 'IGRF coefficients updated 01/01/20',
              'systems': ["GEI", "GEO", "GSE", "GSM", "SM", "MAG"] # More are available, but not listed here.
            },
            {
              'name': 'spacepy',
              'version': spacepy.__version__,
              'version_info': None,
              'systems': systems_spacepy
            },
            {
              'name': 'spacepy-irbem',
              'version': spacepy.__version__,
              'version_info': None,
              'systems': systems_spacepy
            },
            {
              'name': 'spiceypy1',
              'version': spiceypy.__version__,
              'version_info': None,
              'systems': ["GEI", "GEI_TOD", "GEI_MOD", "MEAN_ECLIP", "GEO", "GSE", "GSM", "MAG", "SM"],
              'kernels': ['naif0012.tls', 'rbsp_general011.tf', 'de440s.bsp', 'pck00011.tpc']
            },
            {
              'name': 'spiceypy2',
              'version': spiceypy.__version__,
              'version_info': None,
              'systems': ["GEI", "GEI_TOD", "GEI_MOD", "MEAN_ECLIP", "GEO", "GSE", "GSE_Z_PRIMARY", "GSM", "MAG", "SM"],
              'kernels': ['naif0012.tls', 'rbsp_general012_mod.tf', 'de440s.bsp', 'pck00011.tpc']
            },
            {
              'name': 'sunpy',
              'version': sunpy.__version__,
              'version_info': None,
              'systems': ["GEI", "GEO", "GSE", "GSM", "SM", "MAG"],
              'system_aliases': {
                "GEI": "geocentricearthequatorial",
                "GSE": "geocentricsolarecliptic",
                "GSM": "geocentricsolarmagnetospheric",
                "GEO": "itrs",
                "SM": "solarmagnetic",
                "MAG": "geomagnetic"
              }
            },
            {
              'name': 'sscweb',
              'version': None, # TODO: Use current date?
              'version_info': None,
              'notes': 'https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml',
              'systems': ["GEI", "GEO", "GM", "GSE", "GSM", "J2000", "SM"]
            },
            {
            'name': 'pyspedas',
            'version': None,
            'version_info': None,
            'systems': ["GSE", "GSM", "GEI", "SM", "GEO", "J2000"]
            }
          ]

  if info == True:
    return knowns
  else:
    return [x['name'] for x in knowns]
