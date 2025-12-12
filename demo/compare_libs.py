import numpy as np
import hxform

v = np.array([1., 1., 1.])/np.sqrt(3)
t = [2013, 12, 30, 0, 0, 0]

kwargs = {
  "frame_in": "GEO",
  "frame_out": "GEI",
  "libs_exclude": None
}

results = hxform.compare(v, t, **kwargs)
with open(f'./compare_libs/{kwargs["frame_in"]}2{kwargs["frame_out"]}.txt', 'w') as f:
  f.write(results["log"])
