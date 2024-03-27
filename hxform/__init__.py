
# SpacePy
import warnings
warnings.filterwarnings("ignore", message="Leapseconds.*")

def xprint(msg: str):
  """Print to console and log file based on calling script name.

  The first call from myscript.py will create (or overwrite) myscript.log.

  Subsequent calls will append to myscript.log.

  Example
  --------
  >>> from hxform import xprint
  >>> xprint("Log message")
  """

  import os
  from inspect import stack

  fname = stack()[1][1]
  logfile = os.path.realpath(fname[0:-2]) + "log"

  if not 'counter' in xprint.__dict__:
    xprint.counter = {fname: 0}

  if xprint.counter[fname] == 0:
    if os.path.isfile(logfile):
      os.remove(logfile)

  xprint.counter[fname] += 1
  print(msg)
  with open(logfile, "a") as f:
    f.write(str(msg) + "\n")
