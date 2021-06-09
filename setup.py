from setuptools import setup, find_packages

install_requires = ["numpy"]

setup(
    name='hxform',
    version='0.0.1',
    author='Angel Gutarra-Leon, Bob Weigel, Gary Quaresima',
    author_email='rweigel@gmu.edu',
    packages=find_packages(),
    description='Geophysical coordinate transformations using various libraries',
    install_requires=install_requires
)
