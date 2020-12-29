"""Sheath Model 2D. Main program."""

import numpy as np
from math import sin
import matplotlib.pyplot as plt

from packages.Model.Common.Multi_Particle import MULTI_PARTICLE
from packages.Constants import (PI, AMU, UNIT_CHARGE, EV2J, J2EV)

domain = (0.0, 0.01, -0.025, 0.025)
ptcl = PARTICLE('Ar+')
ptcl.init_ptcl('Ion', 40, +1)

dt = 1e-8
Vdc=100
Vrf=10
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

mp = MULTI_PARTICLE('IEDF')
num_per_launch = 100
Arp = {'name':'Ar+','type':'Ion', 'mass':40.0, 'charge':1.0}
posn = np.random.rand(num_per_launch, 3)
vel = np.random.rand(num_per_launch, 3)
mp.gen_particles(num=num_per_launch, prop=Arp, posn=posn, vel=vel)

cllct = list()
step = 0
for i in range(num_ptcl/num_per_launch):
    # init
    posn = np.zeros((num_per_launch, 3))
    vel = np.zeros((num_per_launch, 3))
    mp.update_posn(posn)
    mp.update_vel(vel)
    phi0 = np.random.rand(num_per_launch)*2.0*PI
    t = 0
    dt = 1e-8
    
    while ptcl.isAlive.any():
        
        # move
        E = Edc + Erf*sin(2*PI*freq*t + phi0)
        t += dt
        mp.posn += mp.vel*dt*mp.isAlive
        mp.vel[:, 0] += E*(mp.charge*UNIT_CHARGE)/(mp.mass*AMU)*dt*mp.isAlive
        step += 1
        if step > max_step:
            break
    
    cllct += mp.vel.T.tolist()[0]
    

fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                       constrained_layout=True)
ax.hist(erg, 2*int(Vrf), density=True)
ax.set_xlim(40, 160)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
fig.savefig(f'Vdc{Vdc}_Vrf{Vrf}.png', dpi=600)
plt.close()
