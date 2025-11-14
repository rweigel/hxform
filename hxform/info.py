def libs():
  return list(lib_info().keys())

def frames(lib):
  """Return list of frames supported by a given library."""
  return lib_info(lib)["frames"]

def lib_info(lib=None):
  """Return list of known libraries.

    lib_info() returns a dictionary with keys of library names and values
    of information about each library.

    lib_info(lib) returns a list of strings of library names. Equivalent to 
    list(lib_info().keys()).
  """
  import os

  import sunpy
  import spacepy
  import spiceypy

  # https://spacepy.github.io/coordinates.html
  frames_spacepy = [
    "ECI2000",
    "ECIMOD",
    "ECITOD",
    "GEI",
    "GSE",
    "GSM",
    "GEO",
    "SM",
    "MAG"
  ]

  kernel_dir = os.path.join(os.path.dirname(__file__), '..', 'kernels')
  kernel_dir = os.path.abspath(kernel_dir)

  infos = {
    'cxform': {
      'name': 'cxform',
      'version': None, # TODO: Use commit hash.
      'version_info': None,
      'frames': [
        "J2000",
        "GEI",
        "GEO",
        "GSE",
        "GSM",
        "SM",
        "MAG",
        "RTN",
        "GSEQ",
        "HEE",
        "HAE",
        "HEEQ"
      ]
    },
    'geopack_08_dp': {
      'name': 'geopack_08_dp',
      'version': None,
      'version_info': 'IGRF coefficients updated 01/01/20',
      # More frames are available, but not listed here.
      'frames': [
        "GEI",
        "GEO",
        "GSE",
        "GSM",
        "SM",
        "MAG"
      ]
    },
    'spacepy': {
      'name': 'spacepy',
      'version': spacepy.__version__,
      'version_info': None,
      'frames': frames_spacepy
    },
    'spacepy-irbem': {
      'name': 'spacepy-irbem',
      'version': spacepy.__version__,
      'version_info': None,
      'frames': frames_spacepy
    },
    'spiceypy1': {
      'name': 'spiceypy1',
      'version': spiceypy.__version__,
      'version_info': None,
      'frames': [
        "GEI",
        "GEI_TOD",
        "GEI_MOD",
        "MEAN_ECLIP",
        "GEO",
        "GSE",
        "GSM",
        "MAG",
        "SM"
      ],
      'kernel': {
        'dir': kernel_dir,
        'files': ['naif0012.tls', 'rbsp_general011.tf', 'de440s.bsp', 'pck00011.tpc']
      }
    },
    'spiceypy2': {
      'name': 'spiceypy2',
      'version': spiceypy.__version__,
      'version_info': None,
      'frames': [
        "GEI",
        "GEI_TOD",
        "GEI_MOD",
        "MEAN_ECLIP",
        "GEO",
        "GSE",
        "GSE_Z_PRIMARY",
        "GSM",
        "MAG",
        "SM"
      ],
      'kernel': {
        'dir': kernel_dir,
        'files': ['naif0012.tls', 'rbsp_general012_mod.tf', 'de440s.bsp', 'pck00011.tpc']
      }
    },
    'sunpy': {
      'name': 'sunpy',
      'version': sunpy.__version__,
      'version_info': None,
      'frames': [
        "GEI",
        "GEO",
        "GSE",
        "GSM",
        "SM",
        "MAG"
      ],
      'system_aliases': {
        "GEI": "geocentricearthequatorial",
        "GSE": "geocentricsolarecliptic",
        "GSM": "geocentricsolarmagnetospheric",
        "GEO": "itrs",
        "SM": "solarmagnetic",
        "MAG": "geomagnetic"
      }
    },
    'sscweb': {
      'name': 'sscweb',
      'version': None, # TODO: Use current date?
      'version_info': None,
      'notes': 'https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml',
      'frames': [
        "GEI",
        "GEO",
        "MAG",
        "GSE",
        "GSM",
        "J2000",
        "SM"
      ],
      'system_aliases': {
        'MAG': 'GM'
      }
    },
    'pyspedas': {
      'name': 'pyspedas',
      'version': None,
      'version_info': None,
      'frames': [
        "GSE",
        "GSM",
        "GEI",
        "SM",
        "GEO",
        "J2000"
      ]
    }
  }

  if lib is not None:
    return infos[lib]

  return infos
