"""
2D PLAsma React Module.

REACT2D contains:
    
"""

import numpy as np
from copy import deepcopy

class REACT2D(object):
    """Define the base tranport module/object."""
    
    def __init__(self, name='Rct2d'):
        """
        Init REACT2D.
        
        name: str, var, name of the REACT2D.
        """
        self.name = name
        
    def from_PLASMA(self, PLA):
        """Init REACT2D."""
        self.Se = np.zeros_like(PLA.ne)
        self.Si = np.zeros_like(PLA.ne)
    
    def to_PLASMA(self, PLA):
        """Copy var to PLASMA2D."""
        PLA.Se = deepcopy(self.Se)
        PLA.Si = deepcopy(self.Si)
    
    def calc_src(self, PLA, ke=2.34e-14):
        """Calc src due to ionization."""
        self.Se = ke * np.power(PLA.Te, 0.59)
        self.Se *= np.exp(-17.8/PLA.Te)
        self.Se *= np.multiply(PLA.ne, PLA.nn)
        self.Si = deepcopy(self.Se)
        self._set_nonPlasma(PLA)

    def _set_nonPlasma(self, PLA):
        """Impose fixed Te on the non-PLAsma materials."""
        for idx, mat in np.ndenumerate(PLA.mesh.mat):
            if mat:
                self.Se[idx] = 0.0
                self.Si[idx] = 0.0