"""Sheath Model 2D. Main program."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from packages.Sheath2D.Sheath2D_ptcl import PARTICLE

domain = (0.0, 0.01, -0.025, 0.025)
ptcl = PARTICLE('Ar+')
ptcl.init_ptcl('Ion', 40, +1)

dt = 1e-6
EF = (0.0, -1.0)
EF = np.asarray(EF)
w_loc = 0.001

num_ptcl = 100
max_step = 100

for i in range(num_ptcl):
    # init
    ptcl.init_posn(domain)
    ptcl.init_uvec(['Zero'])
    ptcl.init_erg(['Zero'])
    print(ptcl.step, ptcl.ang, ptcl.uvec, ptcl.posn, ptcl.speed)
    while ptcl.isAlive:
        # move
        ptcl.move_ptcl(dt, EF)
        print(ptcl.step, ptcl.ang, ptcl.uvec, ptcl.posn, ptcl.speed)
        ptcl.check_bndy(domain, 'Periodic')
        ptcl.check_wafer(w_loc)
        if ptcl.step > max_step:
            ptcl.isAlive=False

