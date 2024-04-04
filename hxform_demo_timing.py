import time

import numpy as np
from numpy import matlib

import hxform
from hxform import xprint # print to console and timing.log

N = 10000

libs = hxform.info.known_libs(info=False)
libs.remove('sscweb')

p_in = np.random.random((N, 3))

t = np.array([2001, 1, 1, 0, 0, 0])
kwargs = {
  'csys_in': 'GSM',
  'csys_out': 'GSE',
  'ctype_in': 'car',
  'ctype_out': 'car'
}

# This is done to avoid the overhead of internal imports
for lib in libs:
  kwargs['lib'] = lib
  p_out = hxform.transform(p_in, t, **kwargs)

xprint(60*'-')
xprint(f'# 1 time value; {N} vectors')
xprint('time      lib')
for lib in libs:
  kwargs['lib'] = lib
  p_out = hxform.transform(p_in, t, **kwargs)
  xprint('{0:2.5f}   {1:s}'.format(hxform.transform.execution_time, lib))
xprint(60*'-')


xprint(60*'-')
xprint(f'# {N} identical time values; {N} different vectors')
xprint('time      lib')
t = matlib.repmat(t, N, 1)
for lib in libs:
  if N > 100 and lib == 'sunpy':
    xprint('skipped   {0:s}'.format(lib))
    continue
  kwargs['lib'] = lib
  p_out = hxform.transform(p_in, t, **kwargs)
  xprint('{0:2.5f}   {1:s}'.format(hxform.transform.execution_time, lib))
xprint(60*'-')


xprint(60*'-')
xprint(f'# {N} different time values; {N} different vectors')
xprint('time      lib')
from datetime import datetime, timedelta
delta = timedelta(hours=1)
t = []
dto = datetime(2001, 1, 1, 0, 0, 0)
for i in range(N):
  dt = dto + i*delta
  t.append([dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second])
t = np.array(t)
for lib in libs:
  if N > 100 and lib == 'sunpy':
    xprint('skipped   {0:s}'.format(lib))
    continue
  if N > 1000 and lib == 'spacepy':
    xprint('skipped   {0:s}'.format(lib))
    continue
  kwargs['lib'] = lib
  p_out = hxform.transform(p_in, t, **kwargs)
  xprint('{0:2.5f}   {1:s}'.format(hxform.transform.execution_time, lib))
xprint(60*'-')

import os
os.rename('hxform_demo_timing.log', f'hxform_demo_timing_N-{N}.log')