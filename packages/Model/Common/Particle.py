"""
Particle.py serves as a data center/hub,
    is shered by all particle-based models and Multi_Particle.py,
    supports one single particle.
"""

import numpy as np

from packages.Constants import (PI, AMU, UNIT_CHARGE, EV2J, J2EV)

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
        posn: arr(3) of float, position of particle, unit in m
        vel: arr(3) of float, velocity of particle, unit in m/s
        erg: float, energy of particle, unit in eV
        """
        self.name = name
        self.ptype = list()
        self.mass = np.array(list())
        self.charge = np.array(list())
        self.isAlive = np.array(list())
        self.posn = np.array(list())
        self.vel = np.array(list())
        self.erg = np.array(list())
     
    def read_species(self, fname):
        """Read in species info from a database (Species.csv)."""
        pass
    
    def select_ptcl(self, sp_name):
        """Select particle from the database and assign the info."""
        pass
    
    def customize_ptcl(self, pname, ptype, mass, charge, isAlive=True):
        """Customize a particle."""
        self.pname = pname
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
    
    def vel2erg(self):
        """Convert velocity to energy."""
        temp = np.power(self.vel, 2)
        temp = np.sum(temp)
        self.erg = 0.5*(self.mass*AMU)*temp*J2EV