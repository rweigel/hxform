from setuptools import dist, setup, find_packages

try:
  # https://stackoverflow.com/a/60740731/18433855
  # This script depends on numpy, which may not be installed.
  # The following line installs it. (try/except because
  # if already installed, an error is thrown.)
  dist.Distribution().fetch_build_eggs(['numpy>=1.23.0'])
except:
  pass

from numpy.distutils.core import setup, Extension

install_requires = ["numpy"]

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
    version='0.0.4',
    author='Angel Gutarra-Leon, Bob Weigel, Gary Quaresima',
    author_email='rweigel@gmu.edu',
    packages=find_packages(),
    description='Heliophysical coordinate transformations using various libraries',
    setup_requires=['numpy'],
    install_requires=install_requires,
    ext_modules=[ext1, ext2]
)
