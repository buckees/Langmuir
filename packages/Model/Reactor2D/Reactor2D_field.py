"""
2D Plasma Field Module.

FIELD2D contains:
    Maxwell Field equation
    Output: Ey
"""

import numpy as np
from math import sqrt
from copy import deepcopy

class FIELD2D(object):
    """Define the Eon Energy Module."""
    
    def __init__(self, name='Field2d'):
        """
        Init the FIELD2D.
        
        name: str, var, name of the FIELD2D.
        """
        self.name = name
    
    def from_PLASMA(self, PLA):
        """Copy var from PLASMA2D."""
        self.pwr_in_tot = deepcopy(PLA.pwr_in_tot)
        
    def to_PLASMA(self, PLA):
        """Copy var to PLASMA2D."""
        PLA.Ey = deepcopy(self.Ey)
    
    def readin_Ey(self,fname):
        """Read in E-Field from external file."""
        Ey = np.fromfile(fname, dtype='f4')
        Ey = np.reshape(Ey, (41, 51))
        Ey = np.flip(Ey, 0)
        self.Ey = Ey
    
    def create_Ey(self, MESH):
        """
        Create Ey for play.
        
        MESH: obj MESH2D()
        """
        self.Ey = MESH.z/MESH.z.sum()*1.0e4
        
    def add_Efunc(self, Efunc, MESH):
        """
        Add Efunc for play.
        
        Efunc: function of space.
        """
        self.Ey = Efunc(MESH.x, MESH.z)
    
    def adjust_E(self, power):
        """Adjust E-field to target desired input power."""
        if self.pwr_in_tot:
            fac = power/self.pwr_in_tot
        else:
            fac = 1.0
        fac = min(1.1, fac)
        self.Ey = sqrt(fac)*self.Ey
        