"""
Field.py serves as a data center/hub,
    controls input and output of E and B,
    stores E and B,
    is shered by all particle-based models,
    supports single and multi particle.
"""

import numpy as np
from math import sin

from packages.Model.Common.Field import FIELD
from packages.Constants import PI
        
class FIELD_SHEATH(FIELD):
    """Create FIELD_SHEATH() object just for sheath model."""
    
    def add_para(self, d_sh, Vdc, Vrf, freq, phi):
        """
        Add field func parameters.
        
        d_sh: float, unit in m, sheath thickness
        Vdc, Vrf: float, unit in V, voltage of DC and RF component
        freq: float, unit in Hz, freq of the sheath voltage
        phi: float, unit in rad, phase
        """
        self.d_sh = d_sh
        self.Vdc = Vdc
        self.Vrf = Vrf
        self.freq = freq
        self.phi = phi
    
    def update_phi(self, phi):
        """Update phi."""
        self.phi = phi
    
    def calc_E(self, t):
        """
        Add field func parameters.
        
        t: float, unit in s, time
        """
        E = np.zeros(3)
        E[1] = -self.Vdc/self.d_sh
        E[1] -= self.Vrf/self.d_sh*sin(2*PI*self.freq*t + self.phi)
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