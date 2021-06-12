import numpy as np
from hxform import hxform as hx

from hxform.xprint import Xprint as Xp
xp = Xp() # Print to console and log file

time1 = [1997,1,1]
time2 = [time1, time1]
input1 = [1., 0., 0.]
input2 = [input1,input1]

if False:
    # time and input can be np.arrays
    # Execution time will be faster in this case.
    time1 = np.array(time1)
    time2 = np.array(time2)
    input1 = np.array(input1)
    input2 = np.array(input2)

# Single time, single vector
output = hx.GSMtoGSE(input1, time1)
xp.xprint(output)

# Multiple times, single vector
output = hx.GSMtoGSE(input1, time2)
xp.xprint(output)

# Single time, multiple vectors
output = hx.GSMtoGSE(input2, time1)
xp.xprint(output)

# Multiple times, multiple vectors
output = hx.GSMtoGSE(input2, time2)
xp.xprint(output)
