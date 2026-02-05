import os
import datetime

import numpy

import hxform

omit = ['sscweb']

def header(msg):
  hxform.xprint(60*'-')
  hxform.xprint(msg)
  hxform.xprint(60*'-')
  hxform.xprint('time      lib')


def run(N):

  t = numpy.array([2001, 1, 1, 0, 0, 0])
  v = numpy.random.random((N, 3))

  # This is done to avoid the overhead of imports inside function (which are)
  # cached after the first call.
  for lib in libs:
    kwargs['lib'] = lib
    vp = hxform.transform(v, t, **kwargs)

  header(f'1 time value; {N} different vectors')
  for lib in libs:
    kwargs['lib'] = lib
    vp = hxform.transform(v, t, **kwargs)
    hxform.xprint('{0:2.5f}   {1:s}'.format(hxform.transform.execution_time, lib))


  header(f'{N} identical time values; {N} different vectors')
  t = numpy.tile(t, (N, 1))
  for lib in libs:
    if N > 100 and lib == 'sunpy':
      hxform.xprint('          {0:s} (skipped b/c N > 100)'.format(lib))
      continue
    kwargs['lib'] = lib
    vp = hxform.transform(v, t, **kwargs)
    hxform.xprint('{0:2.5f}   {1:s}'.format(hxform.transform.execution_time, lib))


  header(f'# {N} different time values; {N} different vectors')
  delta = datetime.timedelta(hours=1)
  t = []
  dto = datetime.datetime(2001, 1, 1, 0, 0, 0)
  for i in range(N):
    dt = dto + i*delta
    t.append([dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second])

  t = numpy.array(t)
  for lib in libs:
    if N > 100 and lib == 'sunpy':
      hxform.xprint('          {0:s} skipped b/c N > 100'.format(lib))
      continue
    if N > 1000 and lib == 'spacepy':
      hxform.xprint('          {0:s} skipped b/c N > 1000'.format(lib))
      continue
    kwargs['lib'] = lib
    vp = hxform.transform(v, t, **kwargs)
    hxform.xprint('{0:2.5f}   {1:s}'.format(hxform.transform.execution_time, lib))

  os.rename('timing.log', f'timing/timing_N-{N}.log')


libs = hxform.info.libs()
libs = list(set(libs) - set(omit))


kwargs = {
  'frame_in': 'GSM',
  'frame_out': 'GSE',
  'ctype_in': 'car',
  'ctype_out': 'car'
}

for N in [100, 1000, 10000]:
  run(N)