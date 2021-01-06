"""
2D Plasma Field Module.

FIELD2D contains:
    Maxwell Field equation
    Output: Ey
"""

import numpy as np
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
        self.x = deepcopy(PLA.mesh.x)
        self.z = deepcopy(PLA.mesh.z)
        
    def to_PLASMA(self, PLA):
        """Copy var to PLASMA2D."""
        PLA.Ey = deepcopy(self.Ey)
    
    def readin_Ey(self,fname):
        """Read in E-Field from external file."""
        Ey = np.fromfile(fname, dtype='f4')
        Ey = np.reshape(Ey, (41, 51))
        Ey = np.flip(Ey, 0)
        self.Ey = Ey
    
    def create_Ey(self):
        """Create Ey for play."""
        self.Ey = self.z/self.z.sum()*1.0e4
        