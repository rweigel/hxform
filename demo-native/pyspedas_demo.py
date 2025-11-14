import os
import numpy as np
os.environ["PYSPEDAS_LOGGING_LEVEL"] = "error"

from pyspedas import cotrans
from pyspedas import time_double

#data_in = np.array([np.array([1., 1., 1.])/np.sqrt(3)])
data_in = np.array([[3.46410162, 3.46410162, 3.46410162]])
#time_in = [time_double('2015-12-30T00:00:00')]
time_in = [time_double('2010-01-01T00:00:00')]
print(time_in)
res = cotrans(time_in=time_in, data_in=data_in, coord_in='GEO', coord_out='GSE')
print(res)