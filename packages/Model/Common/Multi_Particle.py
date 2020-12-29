"""
Multi_Particle.py serves as a data center/hub,
    is shered by all particle-based models,
    uses DataFrame from pandas,
    supports an array of particles.
"""

import numpy as np
import pandas as pd

class MULTI_PARTICLE(object):
    """Create MULTI_PARTICLE() object."""

    def __init__(self, name='Multi_Particle'):
        """
        Init the MULTI_PARTICLE().
        
        name: str, var, name of the MULTI_PARTICLE.
        """
        self.name = name
        
    def add_PARTICLE(self, PARTICLE):
        """Add PARTICLE() object to MULTI_PARTICLE()."""
        pass
    
    def to_array(self, weight):
        """Convert data to Numpy array type."""
        pass
    
    def reinit_posn(self):
        """Re-init positions for all particles."""
        pass
    
    def reinit_vel(self):
        """Re-init velocities for all particles."""
        pass
    
    def to_DataFrame(self):
        """Convert data to DataFrame type."""
        pass
    