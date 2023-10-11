import numpy as np
from hxform import hxform as hx

from hxform.xprint import Xprint as Xp
xp = Xp() # Print to console and log file

time1 = [1997, 1, 1]
#time1 = [2010, 1, 10]
#time1 = [2020, 8, 10, 12, 0, 0]
time1 = [2010, 12, 30, 0, 0, 0]
time2 = [time1, time1]
input1 = [1., 1., 1.]
#input1 = [0.50000000, 0.50000000, 0.50000000]
input2 = [input1, input1]

if False:
    # time and input can be np.arrays
    # Execution time will be faster in this case.
    time1 = np.array(time1)
    time2 = np.array(time2)
    input1 = np.array(input1)
    input2 = np.array(input2)

libs = ['cxform','geopack_08_dp','spacepy','spacepy-irbem']
libs2 = ['cxform','geopack_08_dp','spacepy','spacepy-irbem']
initial = 'GSM'
final = 'GSE'
#initial = 'GEO'
#final = 'GSE'

xp.xprint("Time: {}".format(time1))
xp.xprint(f"Transform: {initial} => {final}")
xp.xprint("Input: {}\n".format(input1))
xp.xprint("Output:\n" + 48*"-")

lengths = []
for lib in libs:
  lengths.append(len(lib))
ml = max(lengths)

outputs = np.full((len(libs),3), fill_value=np.nan)

for i, lib in enumerate(libs):
  #xp.xprint('Using library: {}'.format(lib))
  # Single time, single vector
  output = hx.transform(input1, time1, initial, final, lib=lib)
  outputs[i,:] = np.array(output)
  xp.xprint("{}:{}   {:11.8f} {:11.8f} {:11.8f}".format(lib, (ml-len(lib))*" ",*output))

max, min = np.max(outputs, axis=0), np.min(outputs, axis=0)
xp.xprint("\n")
xp.xprint("max-min:              {:11.8f} {:11.8f} {:11.8f}".format(*(max-min)))
xp.xprint("100*|max-min|/|max|:  {:10.4f}% {:10.4f}% {:10.4f}%".format(*(100*np.abs(max-min)/np.abs(max))))

if False:
    # API Demo
    for lib in libs:
        xp.xprint('Using library: {}'.format(lib))
        # Single time, single vector
        output = hx.GSMtoGSE(input1, time1, lib=lib)
        xp.xprint(output)

        # Multiple times, single vector
        output = hx.GSMtoGSE(input1, time2, lib=lib)
        xp.xprint(output)

        # Single time, multiple vectors
        output = hx.GSMtoGSE(input2, time1, lib=lib)
        xp.xprint(output)

        # Multiple times, multiple vectors
        output = hx.GSMtoGSE(input2, time2, lib=lib)
        xp.xprint(output)
