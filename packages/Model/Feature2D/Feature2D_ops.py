"""Feature Model 2D. Operation parameters."""

class PARAMETER():
    def __init__(self):
        """
        Init the PARTICLE().
        
        name: str, name of the particle.
        """
        self.num_ptcl = None
        self.max_step = None
        self.step_fac = None
        self.max_rflct = None
        self.surf_norm_range = None
        self.surf_norm_mode = None
        self.num_plot = None
        self.prob_rflct = None