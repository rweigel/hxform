# Demo of calling Space directly instead of through hxform.

import spacepy.coordinates as sc
from spacepy.time import Ticktock

from hxform.xprint import Xprint as Xp
xp = Xp() # Print to console and log file

input = sc.Coords([1,0,0], 'GSM', 'car', use_irbem=True)

input.ticks = Ticktock(['1997-01-01T00:00:00'], 'ISO')
xp.xprint(input)

output = input.convert('GSE', 'car')
xp.xprint(output)

input = sc.Coords([1,0,0], 'GSM', 'car', use_irbem=False)

input.ticks = Ticktock(['1997-01-01T00:00:00'], 'ISO')
xp.xprint(input)

output = input.convert('GSE', 'car')
xp.xprint(output)
