import spacepy.coordinates as sc
from spacepy.time import Ticktock

def xprint(msg):
    # Print to console and logfile
    import os
    print(msg);
    logfile = os.path.realpath(__file__)[0:-2] + "log"
    if not os.path.exists(logfile):
        with open(logfile, "w") as f: pass
    with open(logfile, "a") as f: f.write(str(msg) + "\n")

input = sc.Coords([1,0,0], 'GSM', 'car')

input.ticks = Ticktock(['1997-01-01T00:00:00'], 'ISO')
xprint(input)

output = input.convert('GSE', 'car')
xprint(output)

