"""
Multi_Particle.py serves as a data center/hub,
    is shered by all particle-based models,
    uses DataFrame from pandas,
    supports an array of particles.
"""

import numpy as np
import pandas as pd
from math import cos, sin, sqrt, acos
import matplotlib.pyplot as plt
from scipy.stats import cosine

class MULTI_PARTICLE(object):
    """Create MULTI_PARTICLE() object."""

    def __init__(self, name='Multi_Particle'):
        """
        Init the MULTI_PARTICLE().
        
        name: str, var, name of the MULTI_PARTICLE.
        """
        self.name = name
        
    def add_PARTICLE(self):
        """Add PARTICLE() object to MULTI_PARTICLE()."""
        pass
    
    def to_DataFrame(self):
        """Convert data to DataFrame type."""
        pass
    
    def to_array(self):
        """Convert data to Numpy array type."""
        pass
    
    def reinit_posn(self):
        """Re-init positions for all particles."""
        pass
    
    def reinit_vel(self):
        """Re-init velocities for all particles."""
        pass