"""Sheath Model 2D. Main program."""

import os
import glob

import numpy as np
from math import sin
import matplotlib.pyplot as plt
from copy import deepcopy

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Sheath2D.Sheath2D_main import MAIN
from packages.Model.Sheath2D.Sheath2D_field import FIELD_SHEATH
from packages.Model.Common.Particle_Mover import EULER_MOVE, LEAPFROG
from packages.Constants import PI


def Efunc_single(t):
    d_sh = 0.01
    freq = 14e6
    Vdc = 0
    Vrf = 50
    phi = PI/2
    # return min(-Vdc/d_sh - Vrf/d_sh*sin(2*PI*freq*t), 0.0)
    return -Vdc/d_sh - Vrf/d_sh*sin(2*PI*freq*t + phi)

def Efunc_dual(t):
    d_sh = 0.01
    freq1 = 2e6
    freq2 = 14e6
    Vdc = 100
    Vrf1 = 50
    Vrf2 = 10.0
    # return min(-Vdc/d_sh - Vrf1/d_sh*sin(2*PI*freq1*t) - 
    #            Vrf2/d_sh*sin(2*PI*freq2*t), 0.0)
    return -Vdc/d_sh-Vrf1/d_sh*sin(2*PI*freq1*t)-Vrf2/d_sh*sin(2*PI*freq2*t)

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)
ptcl.update_posn(np.zeros(3))
ptcl.update_vel(np.zeros(3))

ptcl2 = PARTICLE()
ptcl2.customize_ptcl('Ion', 40, 1)
ptcl2.update_posn(np.zeros(3))
ptcl2.update_vel(np.zeros(3))
# init field
field = FIELD_SHEATH('Sheath')
field.add_Efunc(Efunc_single)

t, dt = 0.0, 1e-9
time, traj, Efield = list(), list(), list()
traj2 = list()

for i in range(300):
    time.append(t)
    traj.append(ptcl.posn[1])
    traj2.append(ptcl2.posn[1])
    field.update_E(field.calc_E(t))
    Efield.append(field.E)
    posn_next, vel_next = LEAPFROG(ptcl, field, t, dt)
    ptcl.update_posn(posn_next)
    ptcl.update_vel(vel_next)
    posn_next, vel_next = EULER_MOVE(ptcl2, field, t, dt)
    ptcl2.update_posn(posn_next)
    ptcl2.update_vel(vel_next)
    t += dt


fname = 'LEAPFROG'

for i in glob.glob(fname + '.png'):
    os.remove(i)


fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                           constrained_layout=True)
ax.plot(time, traj, 'b-', label='LEAPFROG')
ax.plot(time, traj2, 'r-', label='EULER_FORWARD')
ax.set_title('f=14MHz, dt=1ns')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Position (m)')
ax.legend()
fig.savefig(fname + '.png', dpi=600)
plt.close()

