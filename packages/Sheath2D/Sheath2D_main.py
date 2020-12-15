"""Sheath Model 2D. Main program."""

import numpy as np
from math import sin
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from packages.Sheath2D.Sheath2D_ptcl import PARTICLE
from packages.Constants import PI

domain = (0.0, 0.01, -0.025, 0.025)
ptcl = PARTICLE('Ar+')
ptcl.init_ptcl('Ion', 40, +1)

dt = 1e-6
EF0 = (0.0, -1e4)
EF0 = np.asarray(EF0)
w_loc = 0.001
freq = 13.56e5

num_ptcl = 10000
max_step = 10000

for i in range(num_ptcl):
    # init
    ptcl.init_posn(domain)
    ptcl.init_uvec(['Zero'])
    ptcl.init_erg(['Zero'])
    t = np.random.uniform(0.0, 2*PI)
    while ptcl.isAlive:
        # move
        EF = EF0 + EF0*sin(2*PI*freq*t)
        ptcl.move_ptcl(dt, EF)
        t += dt
        EF = sin(t)
        ptcl.check_bndy(domain, 'Periodic')
        ptcl.check_wafer(w_loc)
        if ptcl.step > max_step:
            ptcl.isAlive=False

erg = [item[0] for item in ptcl.cllct]
plt.hist(erg, 100)