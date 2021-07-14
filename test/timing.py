import time
import numpy as np
from numpy import matlib
from hxform import hxform as hx

import warnings
# Suppress SpacePy warning
warnings.filterwarnings("ignore", message="Leapseconds.*")

from hxform.xprint import Xprint as Xp
xp = Xp() # Print to console and log file

N = 10000
libs = ['spacepy','spacepy-irbem','cxform','geopack_08_dp']

p_in = np.random.random((N,3))
t = np.array([2001,1,1,0,0,0])

xp.xprint(60*'-')
xp.xprint('t.size = (1, 6); p_in.size = ({0:d}, 3)'.format(N))
xp.xprint('time     lib')
for lib in libs:
    tic = time.time()
    p_out = hx.MAGtoGSM(p_in, t, ctype_in='car', ctype_out='car', lib=lib)
    toc = time.time()
    xp.xprint('{0:.4f}   {1:s}'.format(toc-tic, lib))

xp.xprint(60*'-')
xp.xprint('t.size = ({0:d}, 6); p_in.size = ({0:d}, 3)'.format(N))
xp.xprint('time     lib')

t = matlib.repmat(t, N, 1)

for lib in libs:
    tic = time.time()
    p_out = hx.MAGtoGSM(p_in, t, ctype_in='car', ctype_out='car', lib=lib)
    toc = time.time()
    xp.xprint('{0:.4f}   {1:s}'.format(toc-tic, lib))

if False:
    xp.xprint(60*'-')
    xp.xprint('t.size = ({0:d}, 6); p_in.size = (1, 3)'.format(N))
    xp.xprint('time     lib')
    t = np.matlib.repmat(t, N, 3)
    p_in = p_in[0] # Error
    for lib in libs:
        tic = time.time()
        p_out = hx.MAGtoGSM(p_in, t, ctype_in='car', ctype_out='car', lib=lib)
        toc = time.time()
        xp.xprint('{0:.3f}   {1:s}'.format(toc-tic, lib))

xp.xprint(60*'-')
