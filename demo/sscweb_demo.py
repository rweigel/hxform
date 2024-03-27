import re
import requests

import numpy as np
time = np.array([[2000, 1, 1, 0, 0, 0]])

from hxform import hxform as hx

url = "https://sscweb.gsfc.nasa.gov/cgi-bin/CoordCalculator.cgi?"
url += "epoch=2000 100 00:00:00&x=1&y=0.0&z=0.0&lat=&lon=&r=&action=GEI"

r = requests.get(url)
if r.status_code != 200:
  raise Exception(f"Failed to fetch URL: {url}")

text = r.text

# Split the text into lines
lines = text.split("\n")

# Find the start and end of the table
start = 1 + lines.index("          Lat    Long      X       Y       Z   hh.hhhhh ")
end = -1 + [i for i, item in enumerate(lines) if re.search('^ REGION', item)][0]

# Extract the table lines
table_lines = lines[start:end]

# Parse the table into a list of dictionaries
table = []
result = {}
for line in table_lines:
  parts = line.split()
  result[parts[0]] = {
      "Lat": float(parts[1]),
      "Long": float(parts[2]),
      "X": float(parts[3]),
      "Y": float(parts[4]),
      "Z": float(parts[5]),
  }
  if len(parts) > 6:
      result[parts[0]]["hh.hhhhh"] = float(parts[6])

print(result)