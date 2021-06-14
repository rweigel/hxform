import os
from inspect import stack

class Xprint:
    """Print to console and logfile based on calling script name."""

    def __init__(self):
        fname = stack()[1][1]
        logfile = os.path.realpath(fname[0:-2]) + "log"
        if not os.path.exists(logfile):
            with open(logfile, "w") as f: pass
        self.logfile = logfile

    def xprint(self, msg):
        print(msg);
        with open(self.logfile, "a") as f:
            f.write(str(msg) + "\n")
