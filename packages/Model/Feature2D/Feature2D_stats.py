"""Diagnose the feature model."""

class STATS(object):
    """Mesh object."""

    def __init__(self, name='Chemistry'):
        """
        Init the MESH2D.
        
        name: str, var, name of the CHEMISTRY.
        """
        self.launch = dict()
        self.posn = dict()
        self.vel = dict()
        self.erg = dict()
        self.ang = dict()
        self.rflct = dict()
        self.escape = dict()
        self.etch = dict()
        self.term = dict()