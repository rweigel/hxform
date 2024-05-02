import numpy as np
import hxform

v = np.array([1., 1., 1.])/np.sqrt(3)
t = [2010, 12, 30, 0, 0, 0]
libs_exclude = ['sscweb']
#libs_exclude = []
hxform.compare(v, t, 'GEI', 'GEO', libs='all', libs_exclude=libs_exclude)