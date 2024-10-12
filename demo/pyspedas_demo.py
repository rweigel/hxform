import numpy as np
import os
os.environ["PYSPEDAS_LOGGING_LEVEL"] = "error"

from pyspedas import cotrans

from pyspedas import time_double
data_in = [np.array([1., 1., 1.])/np.sqrt(3)]
data_in = [[3.46410162, 3.46410162, 3.46410162]]
time_in = [time_double('2015-12-30/00:00')]

res = cotrans(time_in=time_in, data_in=data_in, coord_in='GSE', coord_out='GSM')
print(res)
