import numpy as np
import hxform

d = 6.0
v = d*np.array([1., 1., 1.])/np.sqrt(3)
t = [2013, 12, 30, 0, 0, 0]
kwargs = {
  "frame_in": "GSE",
  "frame_out": "GSM",
  "rep_in": "car",
  "rep_out": "car",
#  "libs_exclude": ['sscweb']
}

results = hxform.compare(v, t, **kwargs)
with open(f'compare-{kwargs["frame_in"]}2{kwargs["frame_out"]}.txt', 'w') as f:
  f.write(results["log"])

#print(results)
results['dtheta'] = 0.087
dd = d*results['dtheta']*(np.pi/180.0)*6378.16
print(f"If x, y, z are in R_E, then max-min angle corresponds to {dd:.2f} km")