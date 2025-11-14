# See https://github.com/spacepy/spacepy/issues/534

import sys
import spacepy.coordinates as sc
from spacepy.time import Ticktock

if len(sys.argv) > 1 and sys.argv[1] == "1":
	# Computing SM to GSM first seems to be needed
	# in order for GSM to SM to not return nans.
	input = sc.Coords([1,0,0], 'SM', 'car')

	input.ticks = Ticktock(['1997-01-01T00:00:00'], 'ISO')
	#print(input)

	output = input.convert('GSM', 'car')
	#print(output)

input = sc.Coords([1,0,0], 'GSM', 'car')

input.ticks = Ticktock(['1997-01-01T00:00:00'], 'ISO')
print(input)

output = input.convert('SM', 'car')
print(output)


"""
    python spacepy_demo.py 0
gives
    Coords( [[1, 0, 0]] , 'GSM', 'car')
    Coords( [[nan, nan, nan]] , 'SM', 'car')
"""

"""
    python spacepy_demo.py 1
gives
	Coords( [[1, 0, 0]] , 'GSM', 'car')
	Coords( [[0.9001697081193482, 0.0, -0.4355393169213629]] , 'SM', 'car')
"""

