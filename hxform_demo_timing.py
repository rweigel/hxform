import time

import numpy as np
from numpy import matlib

from hxform import hxform as hx
from hxform import xprint # print to console and timing.log

N = 100

libs = hx.known_libs(info=False)
libs.remove('sscweb')

p_in = np.random.random((N, 3))
t = np.array([2001, 1, 1, 0, 0, 0])
kwargs = {
  'csys_in': 'GSM',
  'csys_out': 'GSE',
  'ctype_in': 'car',
  'ctype_out': 'car'
}


xprint(60*'-')
xprint(f'# 1 time value; {N} vectors')
xprint('time     lib')
for lib in libs:
  tic = time.time()
  kwargs['lib'] = lib
  p_out = hx.transform(p_in, t, **kwargs)
  toc = time.time()
  xprint('{0:2.4f}   {1:s}'.format(toc-tic, lib))


xprint(60*'-')
xprint(f'# {N} time value; {N} vectors')
xprint('time     lib')

t = matlib.repmat(t, N, 1)

for lib in libs:
  tic = time.time()
  kwargs['lib'] = lib
  p_out = hx.transform(p_in, t, **kwargs)
  toc = time.time()
  xprint('{0:2.4f}   {1:s}'.format(toc-tic, lib))

xprint(60*'-')
