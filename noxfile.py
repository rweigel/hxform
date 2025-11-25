"""Usage:
Install nox with PBS support (auto-download of Python versions):
  python -m pip install --upgrade 'nox[pbs]'
Test multiple Python versions:
  nox -s tests
Test using current Python version:
  nox -s tests_current
"""

import nox

python = [
  "3.10",
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
  report = f"--junitxml=test/reports/report-{pyver}.xml"
  session.run("pytest", "--capture=tee-sys", report)
