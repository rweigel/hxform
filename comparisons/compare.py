import numpy as np
import hxform

v = np.array([1., 1., 1.])/np.sqrt(3)
t = [2010, 12, 30, 0, 0, 0]
kwargs = {
  "frame_in": "GEO",
  "frame_out": "GSM",
  "rep_in": "car",
  "rep_out": "car",
  "libs_exclude": [] # ['sscweb']
}

results = hxform.compare(v, t, **kwargs)
with open('compare.txt', 'w') as f:
  f.write(results["log"])

results = hxform.compare(v, t, **kwargs)
with open('compare.txt', 'a') as f:
  f.write("\n")
  f.write(results["log"])
