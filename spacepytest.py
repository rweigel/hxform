import spacepy.coordinates as sc
from spacepy.time import Ticktock
import spacepy
import sys

print("spacepy version: ",spacepy.__version__)
print("python version: ", sys.version)

cvals = sc.Coords([[1,2,4],[1,2,2]], 'GSM', 'car')
cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticks
newcoord = cvals.convert('SM', 'car')
print(newcoord)
