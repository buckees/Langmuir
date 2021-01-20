"""
Field.py serves as a data center/hub,
    controls input and output of E and B,
    stores E and B,
    is shered by all particle-based models,
    supports single and multi particle.
"""

import numpy as np

from packages.Model.Common.Field import FIELD
        
class FIELD_SHEATH(FIELD):
    """Create FIELD_SHEATH() object just for sheath model."""

    def add_Efunc(self, Efunc):
        """
        Add Efunc.
        
        Efunc: function, f(t) is function of time.
        """
        self.Efunc = Efunc
    
    def calc_E(self, t):
        """
        Add field func parameters.
        
        t: float, unit in s, time
        """
        E = np.zeros(3)
        E = self.Efunc(t)
        return E
    
    def read_E(self, fname):
        """
        Read in E-field file.
        
        fname: str, filename for E-field.
        """
        self.E_input = np.genfromtxt(fname, delimiter=',')
    
    def interp_E(self, t):
        """
        Interpolate E-field according to E_input.
        
        t: float, unit in s, time
        """
        E = np.zeros(3)
        period = self.E_input[0].max() - self.E_input[0].min()
        E[1] = np.interp(t, self.E_input[0], self.E_input[1], period=period)
        return E