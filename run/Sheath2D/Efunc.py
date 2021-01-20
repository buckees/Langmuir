"""Efunc.py provides built-in E-field function."""

import numpy as np
from math import sin

from packages.Constants import PI
        
class EFUNC(object):
    """Create EFUNC() object just for sheath model."""

    def load(self, OPER):
        """
        Load parameters from obj OPER.
        
        OPER: obj, contains all parameters.
        """
        if OPER.imode_Efunc == 'Single' or 'Dual':
            self.freq = OPER.freq
            self.Vdc = OPER.Vdc
            self.Vrf = OPER.Vrf
            self.d_sh = OPER.d_sh
            if OPER.imode_Efunc == 'Dual':
                self.freq2 = OPER.freq2
                self.Vdc2 = OPER.Vdc2
                self.Vrf2 = OPER.Vrf2
    
    def sgl_freq(self, t):
        """
        E-field function of time with single frquency.
        
        The E-field is negative, only along z-axis.
        E: arr(3) of float, E-field
        """
        E = np.zeros(3)
        E[1] = -self.Vdc/self.d_sh \
               -self.Vrf/self.d_sh*sin(2*PI*self.freq*t)
        E[1] = min(E[1], 0.0)
        return E
    
    def dual_freq(self, t):
        """
        E-field function of time with single frquency.
        
        The E-field is negative, only along z-axis.
        E: arr(3) of float, E-field
        """
        E = np.zeros(3)
        E[1] = -self.Vdc/self.d_sh \
               -self.Vrf/self.d_sh*sin(2*PI*self.freq*t) \
               -self.Vrf2/self.d_sh*sin(2*PI*self.freq2*t)
        E[1] = min(E[1], 0.0)
        return E