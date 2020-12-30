"""Sheath Model 2D. Main program."""

import numpy as np
from math import sin
import matplotlib.pyplot as plt

from packages.Model.Common.Particle import PARTICLE
from packages.Constants import (PI, AMU, UNIT_CHARGE)

# Physics parameters
w_loc = 0.0 # by default

dsh = 0.01
dt = 1e-8
Vdc=100
Vrf=50
# freq = 13.56e6 # T = 70 ns
freq = 1

# Model parameters
num_ptcl = 10000
max_step = 10000

ptcl = PARTICLE('M+')
ptcl.customize_ptcl('Ion', 40, 1)

erg = list()

for i in range(num_ptcl):
    # init
    ptcl.update_state(True)
    ptcl.update_posn(np.array([0.0, dsh, 0.0]))
    ptcl.update_vel(np.zeros(3))
    phi0 = np.random.uniform(0.0, 2*PI)
    t = 0
    dt = 1e-8
    step = 0
    
    while ptcl.isAlive:
        # make sure not over shoot
        if ptcl.vel[1]:
            dt1 = (w_loc - ptcl.posn[1])/ptcl.vel[1]
            if dt1 < dt:
                dt = dt1*1.001
        # move
        EF = np.zeros(3)
        EF[1] = -Vdc/dsh - Vrf/dsh*sin(2*PI*freq*t + phi0)
        ptcl.update_posn(ptcl.posn + ptcl.vel*dt)
        ptcl.update_vel(ptcl.vel + 
                        EF*(ptcl.charge*UNIT_CHARGE)/(ptcl.mass*AMU)*dt)
        t += dt
        step += 1
        if ptcl.posn[1] < w_loc:
            ptcl.isAlive=False
            erg.append(ptcl.vel2erg())
        if step > max_step:
            ptcl.isAlive=False

# plot results
print(f'{num_ptcl} particles are launched.' 
      + f'\n{len(erg)} particles are collected by the wafer.')


fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                       constrained_layout=True)
ax.hist(erg, 2*int(Vrf), density=True)
ax.set_xlim(40, 160)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
fig.savefig(f'Vdc{Vdc}_Vrf{Vrf}.png', dpi=600)
plt.close()
