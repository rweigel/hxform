import numpy as np
from datetime import datetime

import hxform.geopack_08_dp as geopack_08_dp
import hxform.hxform as hx

pos = [1., 0., 0.]
time = [2003.0, 11.0, 20.0, 7.0, 0.0, 0.0]

v_sp = hx.MAGtoGSM(pos, time, ctype_in='car', ctype_out='car', lib='spacepy')
print(v_sp)

dtime = np.array(time[0:5], dtype=np.int32)
output = np.column_stack(geopack_08_dp.transform(pos[0], pos[1], pos[2], "MAGtoGSM", dtime))
print(output)

