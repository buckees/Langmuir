"""Feature Model 2D. Operation parameters."""

class PARAMETER():
    def __init__(self):
        """
        Init the PARTICLE().
        
        name: str, name of the particle.
        """
        self.num_ptcl = None
        self.max_step = None
        self.d_sh = None
        self.wafer_loc = None
        self.dt = None
        self.Vdc = None
        self.Vrf = None
        self.freq = None