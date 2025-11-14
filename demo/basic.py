import hxform

lib = 'spiceypy1'
final = 'GSM'
initial = 'GSE'
time1 = '2010-12-30T00:00:00Z'
time2 = '2010-12-31T01:00:00Z'
input1 = [1., 1., 1.]
input2 = [1., 0., 1.]


time1 = hxform.timelib.iso2ints(time1)
time2 = hxform.timelib.iso2ints(time2)

# Transform one input at one time
hxform.xprint("")
hxform.xprint(f"Transform: {initial} => {final}")
hxform.xprint(f"Time:      {time1}")
hxform.xprint(f"Input:     {input1}")
hxform.xprint(f"Library:   {input1}")

output = hxform.transform(input1, time1, initial, final, lib=lib)
hxform.xprint(f"Output: {output}\n")

# Transform single input at two different times
hxform.xprint("")
hxform.xprint(f"Transform: {initial} => {final}")
hxform.xprint(f"Time:   {[time1, time2]}")
hxform.xprint(f"Input:  {input1}")
output = hxform.transform(input1, [time1, time2], initial, final, lib=lib)
hxform.xprint(f"Output: {output}\n")

# Transform two inputs at same time
hxform.xprint("")
hxform.xprint(f"Transform: {initial} => {final}")
hxform.xprint(f"Time:   {time1}")
hxform.xprint(f"Input:  {[input1, input2]}")
output = hxform.transform([input1, input2], time1, initial, final, lib=lib)
hxform.xprint(f"Output: {output}\n")

# Transform two inputs, each at different times
time2 = [2010, 12, 31, 0, 0, 0]
hxform.xprint("")
hxform.xprint(f"Transform: {initial} => {final}")
hxform.xprint(f"Time:   {[time1, time2]}")
hxform.xprint(f"Input:  {[input1, input2]}")
output = hxform.transform([input1, input2], [time1, time2], initial, final, lib=lib)
hxform.xprint(f"Output: {output}\n")
