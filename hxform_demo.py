import numpy as np
import hxform.hxform as hx

# Print to console and logfile
import os
logfile = os.path.realpath(__file__)[0:-2] + "log"
if os.path.exists(logfile):
    with open(logfile, "w") as f: pass
def xprint(msg):
    import os; print(msg);
    with open(logfile, "a") as f: f.write(str(msg) + "\n")

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
output = hx.GSMtoSM(input1, time1)
xprint(output)

# Multiple times, single vector
output = hx.GSMtoSM(input1, time2)
xprint(output)

# Single time, multiple vectors
output = hx.GSMtoSM(input2, time1)
xprint(output)

# Multiple times, multiple vectors
output = hx.GSMtoSM(input2, time2)
xprint(output)
