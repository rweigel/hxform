import pathlib
import sys
import subprocess
pyver = f"{sys.version_info.major}.{sys.version_info.minor}"
here = pathlib.Path(__file__).parent
for test_file in here.glob("*_test.py"):
  print(f"Running tests in {test_file} with Python {pyver}")
  subprocess.run([sys.executable, str(test_file)], check=False)
