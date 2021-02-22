"""
Particle.py serves as a data center/hub,
    is shered by all particle-based models and Multi_Particle.py,
    supports one single particle.
"""

import numpy as np
from math import sqrt, acos, degrees
import pandas as pd

from packages.Constants import (AMU, J2EV)

class PARTICLE(object):
    """Create PARTICLE() object."""

    def __init__(self):
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
        self.ptype = None
        self.mass = None
        self.charge = None
        self.isAlive = None
        self.posn = np.zeros(3)
        self.vel = np.zeros(3)
     
    def read_species(self, fsp='Species'):
        """
        Read in species info from a database (Species.csv).
        
        Store species information in DataFrame.
        fsp: str, filename for species database
        df_species: DataFrame, store species information
        """
        self.df_species = pd.read_csv(fsp + '.csv', header=0)
    
    def select_ptcl(self, sp_name):
        """
        Select particle from the database and assign the info.
        
        sp_name: str, name of species.
        """
        if sp_name in self.df_species['Name']:
            row = df_species[df_species['Name'] == sp_name].iloc[0]
            self.name = row['Name']
            self.mass = row['Mass']
            self.charge = row['Charge']
        else:
            return f'\n{sp_name} is not found in the databae, "Species.csv".'
    
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
    
    def add_initPosnFunc(self, xFunc):
        """Add init position function."""
        self.initPosnFunc = xFunc
    
    def init_posn(self):
        """Init position."""
        self.posn = self.initPosnFunc()
    
    def add_initVelFunc(self, vFunc):
        """Add init velocity function."""
        self.initVelFunc = vFunc
    
    def init_vel(self):
        """Init velocity."""
        self.vel = self.initVelFunc()
    
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
        