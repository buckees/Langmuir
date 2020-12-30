"""
Particle.py serves as a data center/hub,
    is shered by all particle-based models and Multi_Particle.py,
    supports one single particle.
"""

import numpy as np

class PARTICLE(object):
    """Create PARTICLE() object."""

    def __init__(self, name='Particle'):
        """
        Init the PARTICLE().
        
        name: str, name of the PARTICLE.
        pname: str, name of particle
        ptype: str, type of particle, ['Eon', 'Ion', 'Neut']
        mass: float, mass of particle, unit in AMU
        charge: float, charge of particle
        isAlive: bool, state of particle
        posn: arr(3) of float, position of particle, unit in m
        vel: arr(3) of float, velocity of particle, unit in m/s
        erg: float, energy of particle, unit in eV
        """
        self.name = name
        self.pname = list()
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
    
    def init_ptcl(self, ptype, mass, charge, isAlive=False):
        """Init particle."""
        self.ptype = ptype  # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass  # unit in AMU
        self.charge = charge  # unit in Unit Charge of Electron
        self.posn = np.zeros(2)
        self.enrg = 0.025  # unit in eV, initial as room temperature
        self.uvec = np.zeros(2)
        self.accl = np.zeros(2)
        self.isAlive = False  # indicator for ptcl alive or dead

    def init_posn(self, posn):
        """Init position."""
        self.posn = posn

    def init_vel(self, vel):
        """Init velocity."""
        self.vel = self.vel