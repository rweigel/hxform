import numpy as np
import hxform.geopack_08_dp as geopack_08_dp
import hxform.hxform as cx

def xprint(msg):
    # Print to console and logfile
    import os
    print(msg);
    logfile = os.path.realpath(__file__)[0:-2] + "log"
    if not os.path.exists(logfile):
        with open(logfile, "w") as f: pass
    with open(logfile, "a") as f: f.write(str(msg) + "\n")

libs=['geopack_08_dp','spacepy']

time = (1997,1,1)
input = [1., 0., 0.]

for lib in libs:
    output = cx.GSMtoSM(input, time, lib=lib)
    xprint(input)
    xprint(output)

    times = [time, time]
    output = cx.GSMtoSM(input, times, lib=lib)
    xprint(input)
    xprint(output)

    inputs = [input,input]
    output = cx.GSMtoSM(input, time, lib=lib)
    xprint(input)
    xprint(output)

    times = [time, time]
    inputs = [input,input]
    output = cx.GSMtoSM(input, time, lib=lib)
    xprint(input)
    xprint(output)