"""
Sheath2D_coll.py contains collision information.
"""

class COLLISION(object):
    """Create COLLISION() object just for sheath model."""
    
    def __init__(self, name='Collision'):
        """
        Init the COLLISION.
        
        name:str, name of the collision.
        """
        self.name = name

    def add_func_CollFreq(self, func_CollFreq):
        """
        Add func of collision frquency.
        
        func_CollFreq: function of velocity, f(velocity) returns collision frequency.
        """
        self.func_CollFreq = func_CollFreq
    
    def add_func_ReinitVel(self, func_ReinitVel):
        """
        Add func of reinit velocity.
        
        func_ReinitVel: function of velocity, f(velocity) returns velocity.
        """
        self.func_ReinitVel = func_ReinitVel