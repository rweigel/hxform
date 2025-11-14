import re
import requests

import numpy as np
time = np.array([[2000, 1, 1, 0, 0, 0]])

from hxform import hxform as hx

#url = "https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi?"
#url += "epoch=2000 100 00:00:00&x=1&y=0.0&z=0.0&lat=&lon=&r=&action=GEI"

url = "https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi?"
url += "epoch=2000 100 00:00:00&x=&y=&z=&lat=90.0&lon=0.0=&r=1.0&action=GEO"

r = requests.get(url)
if r.status_code != 200:
  raise Exception(f"Failed to fetch URL: {url}")

text = r.text

# Split the text into lines
lines = text.split("\n")

start = None
end = None
for idx, line in enumerate(lines):
  if re.search('^ Radial distance', line):
    start = idx
  if re.search('^ REGION', line):
    end = idx - 1
    break

if start is None:
  raise Exception(f"Failed to find start of table in URL: {url}. Returned HTML:\n{text}")

# Extract the table lines
R = float(lines[start].split()[2].strip())
table_lines = lines[start+2:end]

# Parse the table into a list of dictionaries
table = []
result = {}
for line in table_lines:
  parts = line.split()
  result[parts[0]] = {
      "R": R,
      "Lat": float(parts[1]),
      "Long": float(parts[2]),
      "X": float(parts[3]),
      "Y": float(parts[4]),
      "Z": float(parts[5]),
  }
  if len(parts) > 6:
      result[parts[0]]["hh.hhhhh"] = float(parts[6])

import json
print(json.dumps(result, indent=2))

import hxform as hx

print(hx.sph2car(1.0, result["GEO"]["Lat"], result["GEO"]["Long"]))