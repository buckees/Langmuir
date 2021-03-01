"""
Particle.py serves as a data center/hub,
    is shered by all particle-based models and Multi_Particle.py,
    supports one single particle.
"""

import numpy as np
from math import sqrt, sin, cos, acos, degrees
import pandas as pd
from scipy.stats import maxwell

from packages.Constants import (AMU, J2EV, EV2J, K2J)

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
        self.name = None
        self.ptype = None
        self.mass = None
        self.charge = None
        self.isAlive = None
        self.posn = np.zeros(3)
        self.vel = np.zeros(3)
     
    def load_database(self, fsp='Species'):
        """
        Read in species info from a database (Species.csv).
        
        Store species information in DataFrame.
        fsp: str, filename for species database
        df_species: DataFrame, store species information
        """
        self.df_species = pd.read_csv(fsp + '.csv', header=0,
                                      index_col=0)
    
    def select_ptcl(self, sp_name):
        """
        Select particle from the database and assign the info.
        
        sp_name: str, name of species.
        """
        if sp_name in self.df_species.index:
            self.name = sp_name
            self.mass = self.df_species.loc[sp_name]['Mass']
            self.charge = self.df_species.loc[sp_name]['Charge']
            self.ptype = self.df_species.loc[sp_name]['Type']
        else:
            print(f'\n"{sp_name}" is not found in database, "Species.csv".')
    
    def customize_ptcl(self, ptype, mass, charge):
        """Customize a particle."""
        self.ptype = ptype
        self.mass = mass
        self.charge = charge
        
    def update_posn(self, posn):
        """Update position."""
        self.posn = posn
    
    def update_vel(self, vel):
        """Update velocity."""
        self.vel = vel
    
    def update_state(self, state):
        """Update state."""
        self.isAlive = state     
    
    def add_xFunc(self, xFunc):
        """Add init position function."""
        self.xFunc = xFunc
    
    def init_posn(self):
        """Init position to xFunc()."""
        self.posn = self.xFunc()
    
    def add_vFunc(self, vFunc):
        """Add init velocity function."""
        self.vFunc = vFunc
    
    def setVel_vFunc(self):
        """Init velocity to vFunc()."""
        self.vel = self.vFunc()
    
    def setVel_norm(self, Tn=300, v0=np.zeros(3)):
        """
        Set velocity to normal distribution.
        
        Tn: float, var, temperature for neutrals in K.
        """
        v0_x, v0_z, v0_y = v0
        sigma = sqrt(Tn*K2J/(self.mass*AMU))  # a = sqrt(kT/m)
        vx = np.random.normal(v0_x, sigma, 1)[0]
        vz = np.random.normal(v0_z, sigma, 1)[0]
        vy = np.random.normal(v0_y, sigma, 1)[0]
        self.vel = np.array([vx, vz, vy])
        # speed = maxwell.rvs(loc=0.0, scale=a, size=1)
        # self.vel = np.zeros(3)
        # mu, sigma = 0.0, 0.1  # default mean and standard deviation
        # theta = np.random.normal(mu, sigma)
        # self.vel[0], self.vel[1] = sin(theta), -cos(theta)
        # self.vel = speed * self.vel
    
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
        