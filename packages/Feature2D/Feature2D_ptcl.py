"""Feature Model 2D. Partile file."""

import numpy as np
from math import cos, sin, sqrt, acos
import matplotlib.pyplot as plt
from scipy.stats import cosine

class LAUNCH(object):
    """Create LAUNCH action."""
    
    def __init__(self, name='Launch'):
        """
        Init the PARTICLE.
        
        name: str, var, name of the PARTICLE.
        """
        self.name = name
        self.sp_name_set = set()
        self.ptcl_list = list()
    
    def readin_sp_database(self, fname):
        """
        Read in species database.
        
        All launched particles must be selected from the database.
        Customized species are not allowed.
        """
        self.sp_name_set = set(['H', 'Ar+'])
        self.sp_set = set()
        pass
    
    def add_ptcl(self, sp_name, flux):
        """Check and add sp to ptcl_list."""
        if sp_name in self.sp_set:
            # find sp in sp_set
            # self.ptcl_list.append(sp)
            pass
        else:
            return f'Error: {sp_name} is not found in the database.'
    
    def pick_ptcl(self):
        """Pick a ptcl to launch."""
        return None

class PARTICLE(object):
    """Create particle object."""

    def __init__(self, name='Particle'):
        """
        Init the PARTICLE.
        
        name: str, var, name of the PARTICLE.
        """
        self.name = name
      
    def init_ptcl(self, ptype, mass, charge, isAlive=False):
        """Init particle."""
        self.ptype = ptype  # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass  # unit in AMU
        self.charge = charge  # unit in Unit Charge of Electron
        self.posn = np.zeros(2)
        self.enrg = 0.025  # unit in eV, initial as room temperature
        self.uvec = np.zeros(2)
        self.accl = np.zeros(2)
        self.isAlive = False  # indicator for ptcl alive or dead

    def init_posn(self, left, right, top):
        """Initialize the position at the top boundary."""
        init_posx = np.random.uniform(left, right)
        self.posn = np.array([init_posx, top])

    def init_uvec(self, idstrb=None):
        """
        Initialize the velocity direction.

        in (x, -z) (half-down quadrant).
        """
        itype = idstrb[0]
        if itype == 'Uniform2D':
            left, right = idstrb[1]/180.0*np.pi, idstrb[2]/180.0*np.pi
            theta = np.random.uniform(left, right)
        elif itype == 'Uniform3D':
            left, right = idstrb[1]/180.0*np.pi, idstrb[2]/180.0*np.pi
            temp_th = np.random.uniform(left, right)
            temp_phi = np.random.uniform(-np.pi, np.pi)
            temp_vz = cos(temp_th)
            temp_vx = sin(temp_th)*cos(temp_phi)
            # temp_vy = sin(temp_th)*sin(temp_phi)
            temp_v = sqrt(temp_vz**2 + temp_vx**2)
            temp_vz = temp_vz/temp_v
            theta = acos(temp_vz)*np.sign(sin(temp_phi))
        elif itype == 'Normal':
            mu, sigma = 0.0, 0.1  # default mean and standard deviation
            mu, sigma = idstrb[1], idstrb[2]
            theta = np.random.normal(mu, sigma)
        elif itype == 'Cosine':
            scale = idstrb[1]/180.0*np.pi
            theta = cosine.rvs(scale=scale, size=1)
        elif itype == 'Mono':
            theta = idstrb[1]/180.0*np.pi

        self.uvec = np.array([sin(theta), -cos(theta)])

    def init_enrg(self, idstrb='Uniform',
                  enrg_min=1e-2, enrg_max=1e4):
        """Initialize the particle energy."""
        if idstrb == 'Uniform':
            enrg = np.random.uniform(enrg_min, enrg_max)

        elif idstrb == 'Normal':
            pass

        elif idstrb == 'Cosine':
            mu, sigma = 0, 0.1  # mean and standard deviation
            pass

        self.enrg = enrg

    def move_ptcl(self, delta_L):
        """Move each partile in a length of delta_L along its v-vector."""
        self.posn += self.uvec*delta_L

    def bndy_check(self, left, right, top, imode='lost'):
        """
        Check the b.c. for moving ptcl.

        make the ptcl dead if it gets beyond the top bdry
        three modes for vertical bdry are available:
            lost, periodic and reflective
        left bdry is alway at 0.0
        """
        if self.posn[1] >= top:
            self.dead = 1
        elif imode == 'lost':
            if not (left < self.posn[0] < right):
                self.dead = 1
        elif imode == 'periodic':
            self.posn[0] = left + ((self.posn[0] - left) % (right - left))
        elif imode == 'reflective':
            pass

        
    
