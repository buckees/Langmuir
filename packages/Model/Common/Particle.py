"""
Particle.py serves as a data center/hub,
    is shered by all particle-based models and Multi_Particle.py,
    supports one single particle.
"""

import numpy as np
from math import sqrt, acos, degrees

from packages.Constants import (AMU, J2EV)

class PARTICLE(object):
    """Create PARTICLE() object."""

    def __init__(self, name='Particle'):
        """
        Init the PARTICLE().
        
        name: str, name of the particle.
        ptype: str, type of particle, ['Eon', 'Ion', 'Neut']
        mass: float, mass of particle, unit in AMU
        charge: float, charge of particle
        isAlive: bool, state of particle
        posn: arr(3) of float, position in (x, z, y), unit in m
        vel: arr(3) of float, velocity in (x, z, y), unit in m/s
        """
        self.name = name
        self.ptype = list()
        self.mass = np.array(list())
        self.charge = np.array(list())
        self.isAlive = np.array(list())
        self.posn = np.array(list())
        self.vel = np.array(list())
     
    def read_species(self, fname):
        """Read in species info from a database (Species.csv)."""
        pass
    
    def select_ptcl(self, sp_name):
        """Select particle from the database and assign the info."""
        pass
    
    def customize_ptcl(self, ptype, mass, charge, isAlive=True):
        """Customize a particle."""
        self.ptype = ptype
        self.mass = mass
        self.charge = charge
        self.isAlive = isAlive
        self.posn = np.zeros(3)
        self.vel = np.zeros(3)
        
    def update_posn(self, posn):
        """Update position."""
        self.posn = posn
    
    def update_vel(self, vel):
        """Update velocity."""
        self.vel = vel
    
    def update_state(self, state):
        """Update state."""
        self.isAlive = state     
    
    def vel2speed(self):
        """Convert velocity to and return speed and uvec."""
        speed = sqrt(np.sum(self.vel**2))
        uvec = self.vel/speed
        return speed, uvec
    
    def vel2erg(self, nvec=np.array([0, -1, 0])):
        """
        Convert velocity to energy.
        
        nvec: arr(3) of float, normal vector.
        erg: float, energy of the particle.
        theta: float, theta w.r.t. nvec.
        """
        speed, uvec = self.vel2speed()
        erg = 0.5*(self.mass*AMU)*speed**2
        erg *= J2EV
        theta = acos(np.dot(uvec, nvec))*np.sign(uvec)[0]
        theta = degrees(theta)
        return erg, theta
    
    def move_in_space(self, dL):
        """
        Move the particle by a spacial step, dL, regardless of speed.
        
        Assume no field at all.
        dL: float, space step, unit in m.
        """
        speed, uvec = self.vel2speed()
        self.posn += uvec*dL
        
    def move_in_time(self, dt):
        """
        Move the particle by a time step, dt.
        
        Assume no field at all.
        dt: float, time step, unit in s.
        """
        self.posn += self.vel*dt
        