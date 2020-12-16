"""Sheath Model 2D. Main program."""

import numpy as np
from math import sin
import matplotlib.pyplot as plt

from packages.Sheath2D.Sheath2D_ptcl import PARTICLE
from packages.Constants import PI

domain = (0.0, 0.01, -0.025, 0.025)
ptcl = PARTICLE('Ar+')
ptcl.init_ptcl('Ion', 40, +1)

dt = 1e-8
Vdc=100
Vrf=50
Edc = (0.0, -Vdc/0.01)
# Erf = (0.0, -1e3)
Erf = (0.0, -Vrf/0.01)
Edc = np.asarray(Edc)
Erf = np.asarray(Erf)
w_loc = 0.00
# freq = 13.56e6 # T = 70 ns
freq = 1

num_ptcl = 10000
max_step = 10000

for i in range(num_ptcl):
    # init
    ptcl.init_posn(domain)
    ptcl.init_uvec(['Zero'])
    ptcl.init_erg(['Zero'])
    t = np.random.uniform(0.0, 2*PI)
    dt = 1e-8
    
    while ptcl.isAlive:
        
        if ptcl.speed:
            dt1 = abs(ptcl.posn[1] - w_loc)/ptcl.speed
            if dt1 < dt:
                dt = dt1*1.01
        # move
        E = Edc + Erf*sin(2*PI*freq*t)
        E1 = Edc + Erf*sin(2*PI*freq*(t+dt))
        ptcl.move_ptcl(dt, E, E1)
        t += dt
        EF = sin(t)
        ptcl.check_bndy(domain, 'Periodic')
        ptcl.check_wafer(w_loc)
        if ptcl.step > max_step:
            ptcl.isAlive=False


erg = [item[0] for item in ptcl.cllct]
print(len(erg))


fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                       constrained_layout=True)
ax.hist(erg, 2*int(Vrf), density=True)
ax.set_xlim(40, 160)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
fig.savefig(f'Vdc{Vdc}_Vrf{Vrf}.png', dpi=600)
plt.close()
