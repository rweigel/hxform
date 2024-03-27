import spacepy.coordinates as sc
from spacepy.time import Ticktock
from hxform import xprint # print to console and spacepy_demo.log

input = sc.Coords([1,0,0], 'GSM', 'car', use_irbem=True)

input.ticks = Ticktock(['1997-01-01T00:00:00'], 'ISO')
xprint(input)

output = input.convert('GSE', 'car')
xprint(output)
