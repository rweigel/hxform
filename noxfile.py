"""Usage:
Install nox with PBS support (auto-download of Python versions):
  python -m pip install --upgrade 'nox[pbs]'
Test multiple Python versions:
  nox -s tests
Test using current Python version:
  nox -s tests_current
"""

import nox

# Can't use earlier versions b/c of bugs in some libs
python = [
  "3.11",
  "3.12",
  "3.13",
  "3.14"
]

# Run using specified Python versions
@nox.session(python=python)
def tests(session):
  run(session)


# Run using the current interpreter
@nox.session
def tests_current(session):
  run(session)


def run(session):
  """Run the test suite using the current Python interpreter."""
  session.install("pip", "setuptools", "wheel", "pytest")
  session.install("-e", ".")
  pyver = session.name
  report = f"--junitxml=test/log/{pyver}.xml"
  try:
    session.run("pytest", "--capture=tee-sys", report)
  except:
    copy_logs(pyver)
    raise


def copy_logs(pyver):
  import shutil
  import pathlib
  dest_dir = pathlib.Path("test/log")
  dest_dir.mkdir(parents=True, exist_ok=True)
  for src in pathlib.Path("test").glob("*.log"):
    dest = dest_dir / f"{src.stem}.{pyver}.log"
    shutil.move(src, dest)

if __name__ == "__main__":
  copy_logs("10")