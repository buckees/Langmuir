"""Feature Model 2D. Operation parameters."""

class PARAMETER():
    def __init__(self):
        """
        Init the PARTICLE().
        
        name: str, name of the particle.
        """
        self.num_iter = None

    def add_PWRfunc(self, PWR_func):
        """
        Add power function of time.

        PWR_func: Power function of time.
        """
        self.PWRfunc = PWR_func

    def update_pwr(self, t):
        """
        update target power at time, t.

        t: float, unit in s, time.
        """
        self.pwr = self.PWRfunc(t)
