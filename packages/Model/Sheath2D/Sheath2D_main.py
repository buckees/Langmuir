"""Sheath Model 2D. Main program."""

import numpy as np
import matplotlib.pyplot as plt

from packages.Model.Common.Multi_Particle import MULTI_PARTICLE
from packages.Constants import (PI, AMU, UNIT_CHARGE, EV2J, J2EV)

# domain = (0.0, 0.01, -0.025, 0.025)
# ptcl = PARTICLE('Ar+')
# ptcl.init_ptcl('Ion', 40, +1)

dt = 1e-8
sheath_thick = 0.01
Vdc=100
Vrf=0
Edc = -Vdc/sheath_thick
Erf = -Vrf/sheath_thick
w_loc = 0.00
# freq = 13.56e6 # T = 70 ns
freq = 1

num_ptcl = 10000
max_step = 1000

mp = MULTI_PARTICLE('IEDF')
num_per_launch = 1000
Arp = {'name':'Ar+','type':'Ion', 'mass':40.0, 'charge':1.0}
posn = np.random.rand(num_per_launch, 3)
vel = np.random.rand(num_per_launch, 3)
mp.gen_particles(num=num_per_launch, prop=Arp, posn=posn, vel=vel)

cllct = list()
for i in range(int(num_ptcl/num_per_launch)):
    print(f'interation {i+1}')
    # init
    posn = np.zeros((num_per_launch, 3))
    posn[:, 1] = sheath_thick
    vel = np.zeros((num_per_launch, 3))
    mp.update_posn(posn)
    mp.update_vel(vel)
    mp.update_state(np.ones(num_per_launch, dtype=bool))
    
    phi0 = np.random.rand(num_per_launch)*2.0*PI
    t = 0
    dt = 1e-8
    step = 0
    
    while mp.isAlive.any():
        # move
        E = np.zeros((num_per_launch, 3))
        E[:, 1] = Edc + Erf*np.sin(2*PI*freq*t + phi0)
        t += dt
        
        temp_isAlive = mp.isAlive.reshape((mp.isAlive.size, 1))
        mp.posn += mp.vel*dt*temp_isAlive
        temp_chmass = (mp.charge*UNIT_CHARGE)/(mp.mass*AMU)
        temp_chmass = temp_chmass.reshape((temp_chmass.size, 1))
        mp.vel += dt*E*temp_chmass*temp_isAlive
        step += 1
        # check step
        if step > max_step:
            break
        # check wafer
        mp.update_state(np.invert(mp.posn[:, 1] < w_loc))
    
    mp.vel2erg()
    cllct += mp.erg.tolist()
    

fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                       constrained_layout=True)
ax.hist(cllct, density=True)
# ax.set_xlim(40, 160)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
fig.savefig(f'Vdc{Vdc}_Vrf{Vrf}.png', dpi=600)
plt.close()
