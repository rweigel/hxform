import sys
from setuptools import setup, find_packages

# https://numpy.org/doc/stable/reference/distutils_status_migration.html#distutils-status-migration
# Ensure the Python version is compatible
if sys.version_info >= (3, 12, 0):
  raise RuntimeError("Python 3.11 or lower is required due to numpy.distutils deprecation.")

try:
  # https://stackoverflow.com/a/60740731/18433855
  # This script depends on numpy, which may not be installed.
  # The following line installs it. (try/except because
  # if already installed, an error is thrown.)
  # This no longer seems to work. NumPy must be installed before running setup.py.
  #dist.Distribution().fetch_build_eggs(['numpy'])
  import numpy
except:
  raise ImportError("NumPy must be installed before running setup.py")

from numpy.distutils.core import setup, Extension

# sunpy>=7.0.0 due to https://github.com/sunpy/sunpy/pull/8193
# pyspedas>=1.7.28 due to https://github.com/spedas/pyspedas/issues/1207
# spacepy>=0.3.0 due to addition of native python transforms and fix of
# https://github.com/spacepy/spacepy/issues/534; see also
# https://github.com/spacepy/spacepy/pull/536
install_requires = [
  "numpy==1.26.4",
  'sunpy==7.0.0',
  'pyspedas==1.7.28',
  'spacepy==0.6.0',
  'spiceypy==6.0.0',
  'python-dateutil==2.9.0.post0'
]

try:
  # Will work if utilrsw was already installed, for example via pip install -e .
  import utilrsw
except:
  install_requires.append("utilrsw @ git+https://github.com/rweigel/utilrsw")

# https://gist.github.com/johntut/1d8edea1fd0f5f1c057c
# https://github.com/PyCOMPLETE/pypkgexample
# https://numpy.org/devdocs/f2py/distutils.html
ext1 = Extension(
                'hxform.geopack_08_dp',
                sources = [
                            'src/geopack-2008/Geopack-2008_dp_wrapper.for',
                            'src/geopack-2008/Geopack-2008_dp.for',
                            'src/geopack-2008/T96_01.for'
                        ])

ext2 = Extension('hxform.cxform_wrapper',
                sources = [
                            'src/cxform/cxform_wrapper.c',
                            'src/cxform/cxform-manual.c',
                            'src/cxform/cxform-auto.c'
                        ])

setup(
    name='hxform',
    version='0.0.5',
    author='Angel Gutarra-Leon, Bob Weigel, Gary Quaresima',
    author_email='rweigel@gmu.edu',
    packages=find_packages(),
    description='Heliophysical coordinate transformations using various libraries',
    setup_requires=['numpy'],
    install_requires=install_requires,
    ext_modules=[ext1, ext2]
)
