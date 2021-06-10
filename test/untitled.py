import numpy as np
from datetime import datetime
import sys
sys.path.insert(0,'Users/Angel/hxform')
import hxform.geopack_08_dp as geopack_08_dp
import hxform.hxform as hx
sys.path.pop(0)

pos = [1., 0., 0.]
# the error occured with the difference in time formats.
# it'll probably be best to make the conversion in hxform.py 
# space takes [year, month, day, h, m s]
time = [2003.0, 11.0, 20.0, 7.0, 0.0, 0.0]
# geopack takes [year, day of year, h, m, s]
dtime = np.array([2003,324,7,0,0],dtype=np.int32)
# you had dtype = np.array(time[0:5],dtype=np.int32)

v_sp = hx.MAGtoGSM(pos, time, ctype_in='car', ctype_out='car', lib='spacepy')
print(v_sp)

output = np.column_stack(geopack_08_dp.transform(pos[0], pos[1], pos[2], "MAGtoGSM", dtime))
print(output)
