import pathlib
import warnings
from setuptools import setup, Extension, find_packages


# Read version info in hxform/version.py without importing the package
main_ns = {}
ver_path = pathlib.Path(__file__).with_name('hxform') / 'version.py'
with ver_path.open() as vf:
  exec(vf.read(), main_ns)

# Read README.md for long description
readme_path = pathlib.Path(__file__).with_name('README.md')
long_description = readme_path.read_text(encoding='utf-8') if readme_path.exists() else ''

# See issues/issues.md for reason for minimum SunPy, SpacePy, and PySPEDAS versions.
# Note that numpy>=1.26 also appears in pyproject.toml and both should be same.
install_requires = [
  'numpy>=1.26',
  'sunpy>=7.0.3',
  'spacepy>=0.6.0',
  'spiceypy>=6.0.0',
  'pyspedas>=1.7.28',
  'python-dateutil>=2.9.0',
]

try:
  # Will work if utilrsw was already installed, for example via pip install -e .
  import utilrsw
except:
  install_requires.append("utilrsw @ git+https://github.com/rweigel/utilrsw")

cxform_sources = [
    'src/cxform/cxform_wrapper.c',
    'src/cxform/cxform-manual.c',
    'src/cxform/cxform-auto.c',
]
cxform_ext = Extension('hxform.cxform_wrapper',sources=cxform_sources)

def build_geopack():
  """Try to build Geopack using numpy.f2py and copy the produced shared library
  into the `hxform` package directory.

  On failure return error message.
  """
  import os
  import sys
  import glob
  import shutil
  import tempfile

  def compile(src_dir):
    import subprocess

    cmd = [
      sys.executable,
      '-m',
      'numpy.f2py',
      '-c',
      str(src_dir / 'Geopack-2008_dp_wrapper.for'),
      str(src_dir / 'Geopack-2008_dp.for'),
      str(src_dir / 'T96_01.for'),
      '-m', 'geopack_08_dp'
    ]
    print('Executing: ', ' '.join(cmd))
    proc = subprocess.run(cmd, cwd=str(out_dir), capture_output=True, text=True)
    if proc.returncode != 0:
      emsg = f'f2py failed with code {proc.returncode}'
      emsg += f'\nstdout:\n{proc.stdout}'
      emsg += f'\nstderr:\n{proc.stderr}'
      return emsg

    return None

  def copy(src_dir, out_dir):

      # Find the produced extension file (platform-dependent suffix)
      patterns = [str(out_dir / 'geopack_08_dp*.so'),
                  str(out_dir / 'geopack_08_dp*.pyd'),
                  str(out_dir / 'geopack_08_dp*.*')
                ]
      found = []
      for pat in patterns:
        found.extend(glob.glob(pat))

      if not found:
        emsg = f'No Geopack build artifact in list {patterns}'
        raise emsg

      # Copy the first matching artifact into hxform/ as geopack_08_dp<suffix>
      target_dir = pathlib.Path(__file__).with_name('hxform')
      target_dir.mkdir(parents=True, exist_ok=True)
      src_art = pathlib.Path(found[0])
      dest = target_dir / src_art.name
      print('Copying', src_art, '->', dest)
      shutil.copy2(src_art, dest)
      print('Geopack extension copied to', dest)

      return None

  src_dir = pathlib.Path(__file__).with_name('src') / 'geopack-2008'
  if not src_dir.exists():
    raise "Directory src/geopack-2008 does not exist."

  # get terminal width with a sensible fallback
  try:
    term_width = shutil.get_terminal_size(fallback=(80, 24)).columns
  except Exception:
    term_width = int(os.environ.get('COLUMNS', '80')) if os.environ.get('COLUMNS') else 80

  print(3*"\n")
  print(term_width*'-')
  with tempfile.TemporaryDirectory() as out_dir:
    out_dir = pathlib.Path(out_dir)

    emsg = compile(src_dir)
    if emsg:
      return emsg

    try:
      copy(src_dir, out_dir)
    except Exception as e:
      print(term_width*'-')
      return str(e)

  print(term_width*'-')
  print(3*"\n")

  return None

emsg = build_geopack()
if emsg:
  warnings.warn(emsg)

setup(
  name='hxform',
  version=main_ns.get('__version__', '0.0.0'),
  author='Bob Weigel, Angel Gutarra-Leon, and Gary Quaresima',
  author_email='rweigel@gmu.edu',
  packages=find_packages(),
  description='Heliophysical coordinate transformations using various libraries',
  long_description=long_description,
  long_description_content_type='text/markdown',
  include_package_data=True,
  package_data={'': ['README.md']},
  install_requires=install_requires,
  ext_modules=[cxform_ext],
)

