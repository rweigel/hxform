import numpy as np
import hxform
from hxform import xprint

# time and input can be NumPy arrays
time1 = hxform.timelib.iso2ints('2010-12-30T00:00:00Z')
time2 = hxform.timelib.iso2ints('2010-12-31T01:00:00Z')
print(time1, time2)
input1 = [1., 1., 1.]
input2 = [1., 0., 1.]
initial = 'GSE'
final = 'GSM'

lib = 'sscweb'

# Transform one input at one time
xprint("")
xprint(f"Transform: {initial} => {final}")
xprint(f"Time:   {time1}")
xprint(f"Input:  {input1}")

output = hxform.transform(input1, time1, initial, final, lib=lib)
xprint(f"Output: {output}\n")

# Transform single input at two different times
xprint("")
xprint(f"Transform: {initial} => {final}")
xprint(f"Time:   {[time1, time2]}")
xprint(f"Input:  {input1}")
output = hxform.transform(input1, [time1, time2], initial, final, lib=lib)
xprint(f"Output: {output}\n")

# Transform two inputs at same time
xprint("")
xprint(f"Transform: {initial} => {final}")
xprint(f"Time:   {time1}")
xprint(f"Input:  {[input1, input2]}")
output = hxform.transform([input1, input2], time1, initial, final, lib=lib)
xprint(f"Output: {output}\n")

# Transform two inputs, each at different times
time2 = [2010, 12, 31, 0, 0, 0]
xprint("")
xprint(f"Transform: {initial} => {final}")
xprint(f"Time:   {[time1, time2]}")
xprint(f"Input:  {[input1, input2]}")
output = hxform.transform([input1, input2], [time1, time2], initial, final, lib=lib)
xprint(f"Output: {output}\n")
