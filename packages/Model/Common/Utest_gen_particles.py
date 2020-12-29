"""
Test the gen_particles() in the Multi_Particles().
"""

import numpy as np

from Multi_Particle import MULTI_PARTICLE

mp = MULTI_PARTICLE()
num = 100
Arp = {'name':'Ar+','type':'Ion', 'mass':40.0, 'charge':1.0}
posn = np.random.rand(num, 3)
vel = np.random.rand(num, 3)
mp.gen_particles(num=100, prop=Arp, posn=posn, vel=vel)

print(mp)