from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension

install_requires = ["numpy"]

# https://gist.github.com/johntut/1d8edea1fd0f5f1c057c
# https://github.com/PyCOMPLETE/pypkgexample
# https://numpy.org/devdocs/f2py/distutils.html
ext1 = Extension(
                name = 'geopack_08_dp',
                sources = [
							'src/geopack-2008/Geopack-2008_dp_wrapper.for',
							'src/geopack-2008/Geopack-2008_dp.for',
							'src/geopack-2008/T96_01.for'
						])

ext2 = Extension('cxformv',
                sources = [
                            'src/cxform/cxformv.c',
                            'src/cxform/cxform-manual.c',
                            'src/cxform/cxform-auto.c'
                        ])

setup(
    name='hxform',
    version='0.0.2',
    author='Angel Gutarra-Leon, Bob Weigel, Gary Quaresima',
    author_email='rweigel@gmu.edu',
    packages=find_packages(),
    description='Heliophysical coordinate transformations using various libraries',
	setup_requires=['numpy'],
    install_requires=install_requires,
    ext_modules=[ext1, ext2]
)

# The result of f2py compilation is a libarary file that is
# placed in the same directory as setup.py. This copies the 
# library file into the package directory. It seems like there
# should be a way to pass an output directory for the the
# library in either Extension or setup.
import os
import glob
import shutil
for file in glob.glob("*so"):
        if os.path.exists(os.path.join("hxform", file)):
            os.remove(os.path.join("hxform", file))
        shutil.move(file, "hxform")
