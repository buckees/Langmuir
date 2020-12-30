"""
Field.py serves as a data center/hub,
    controls input and output of E and B,
    stores E and B,
    is shered by all particle-based models,
    supports single and multi particle.
"""

import numpy as np

from packages.Constants import (PI, AMU, UNIT_CHARGE, EV2J, J2EV)

class FIELD(object):
    """Create FIELD() object."""

    def __init__(self, name='Field'):
        """
        Init the FIELD().
        
        name: str, name of the field.
        E: arr(3) of float, E-field in (x, z, y), unit in V/m
        B: arr(3) of float, B-field in (x, z, y), unit in T
        """
        self.name = name
        self.E = np.zeros(3)
        self.B = np.zeros(3)
    
    def update_E(self, E):
        """Update E-field."""
        self.E = E
    
    def add_Efunc(self, Efunc):
        """
        Add function for E-field.
        
        Efunc: function which determines E-field.
            Efunc is defined or given by the end user.
            or Efunc can be imported from reactor model.
        """
        self.Efunc = Efunc