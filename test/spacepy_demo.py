import spacepy.coordinates as sc
from spacepy.time import Ticktock

input = sc.Coords([1,0,0], 'GSM', 'car')

input.ticks = Ticktock(['1997-01-01T00:00:00'], 'ISO')
print(input)

output = input.convert('GSE', 'car')
print(output)

