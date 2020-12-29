"""
Multi_Particle.py serves as a data center/hub,
    is shered by all particle-based models,
    uses DataFrame from pandas,
    supports an array of particles.
"""

import numpy as np
import pandas as pd

class MULTI_PARTICLE(object):
    """Create MULTI_PARTICLE() object."""

    def __init__(self, name='Multi_Particle'):
        """
        Init the MULTI_PARTICLE().
        
        name: str, var, name of the MULTI_PARTICLE.
        """
        self.name = name
        
    def add_PARTICLE(self, PARTICLE):
        """Add PARTICLE() object to MULTI_PARTICLE()."""
        pass
    
    def to_array(self, weight):
        """Convert data to Numpy array type."""
        pass
    
    def gen_particles(self, num, prop, posn, vel):
        """
        Generate multi particles of a single species directly.
        
        num: int, num of particles.
        prop: dict, properties of the particle.
        posn: numpy array, (num, 3), position of all particles.
        vel: numpy array, (num, 3), velocity of all particles
        """
        self.num = num
        self.pname = [prop['name'] for i in range(num)]
        self.ptype = [prop['type'] for i in range(num)]
        self.mass = np.ones(num)*prop['mass']
        self.charge = np.ones(num)*prop['charge']
        self.isAlive = np.ones(num, dtype=bool)
        self.posn = posn
        self.vel = vel
        return print(f'{num} of particles are generated!')
    
    def update_posn(self, posn):
        """Update positions for all particles."""
        pass
    
    def update_vel(self, vel):
        """Update velocities for all particles."""
        pass
    
    def to_DataFrame(self):
        """Convert data to DataFrame type."""
        pass

if __name__ == '__main__':
    mp = MULTI_PARTICLE()
    num = 100
    Arp = {'name':'Ar+','type':'Ion', 'mass':40.0, 'charge':1.0}
    posn = np.random.rand(num, 3)
    vel = np.random.rand(num, 3)
    mp.gen_particles(num=100, prop=Arp, posn=posn, vel=vel)