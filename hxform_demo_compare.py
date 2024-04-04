import numpy as np
import hxform
from hxform import hxform as hx
from hxform import xprint

# time and input can be NumPy arrays
time1 = [1997, 1, 1]
time1 = [2010, 12, 30, 0, 0, 0]
input1 = [1., 1., 1.]

libs = hxform.info.known_libs(info=False)

initial = 'GSE'
final = 'GSM'

xprint(f"Time: {time1}")
xprint(f"Transform: {initial} => {final}")
xprint(f"Input: {input1}\n")
xprint("Output:\n" + 48*"-")

ml = max([len(lib) for lib in libs]) # Max length of lib name

outputs = np.full((len(libs), 3), fill_value=np.nan)

# Single time, single vector
for i, lib in enumerate(libs):
  output = hx.transform(input1, time1, initial, final, lib=lib)
  outputs[i,:] = np.array(output)
  xprint("{} {}   {:11.8f} {:11.8f} {:11.8f}".format(lib, (ml-len(lib)+5)*" ",*output))

max, min = np.max(outputs, axis=0), np.min(outputs, axis=0)
xprint("\n")
xprint("max-min:              {:11.8f} {:11.8f} {:11.8f}".format(*(max-min)))
xprint("100*|max-min|/|max|:   {:10.4f}% {:10.4f}% {:10.4f}%".format(*(100*np.abs(max-min)/np.abs(max))))
