"""
Particle.py serves as a data center/hub,
    is shered by all particle-based models and Multi_Particle.py,
    supports one single particle.
"""

import numpy as np
import pandas as pd

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
        else:
            print(f'\n"{sp_name}" is not found in database, "Species.csv".')

ptcl = PARTICLE()
ptcl.load_database()
ptcl.select_ptcl('Ar+')
