"""Feature Model 2D. Partile file."""

import numpy as np
from math import cos, sin, sqrt, acos, copysign, radians, degrees
from scipy.stats import cosine
import matplotlib.pyplot as plt

from packages.Constants import (PI, AMU, UNIT_CHARGE, EV2J, J2EV)

class PARTICLE(object):
    """Create particle object."""

    def __init__(self, name='Particle'):
        """
        Init the PARTICLE.
        
        name: str, var, name of the PARTICLE.
        """
        self.name = name
      
    def init_ptcl(self, ptype, mass, charge, isAlive=True):
        """Init particle."""
        self.ptype = ptype  # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass  # unit in AMU
        self.charge = charge  # unit in Unit Charge of Electron
        self.posn = np.zeros(2)
        self.erg = 0.025  # unit in eV, initial as room temperature
        self.ang = 0 # angle w.r.t. (0, -1)
        self._ang2uvec()
        self.accl = np.zeros(2)
        self.isAlive = True  # indicator for ptcl alive or dead
        self.step = 0  # count the computing step
        self.cllct = list()  # collect the ptcl hitting the wafer
        self._erg2speed()

    def _erg2speed(self):
        """Convert energy (eV) to speed (m/s)."""
        self.speed = sqrt(2.0*(self.erg*EV2J)/(self.mass*AMU))
        self.vel = self.speed*self.uvec
    
    def _speed2erg(self):
        """Convert speed (m/s) to energy (eV)."""
        self.erg = 0.5*(self.mass*AMU)*self.speed**2*J2EV
    
    def _uvec2ang(self):
        """Convert uvec to ang (in degree) w.r.t. (0, -1)."""
        temp = copysign(1.0, self.uvec[0])*acos(self.uvec[1])
        self.ang = degrees(temp)
    
    def _ang2uvec(self):
        """Convert ang (in degree) w.r.t. (0, -1) to uvec."""
        temp = radians(self.ang)
        self.uvec = np.array([sin(temp), -cos(temp)])

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

    def init_erg(self, idstrb=None):
        """Initialize the particle energy."""
        itype = idstrb[0]
        
        if itype == 'Uniform':
            erg_min, erg_max = idstrb[1:]
            erg = np.random.uniform(erg_min, erg_max)

        elif itype == 'Normal':
            pass

        elif itype == 'Cosine':
            mu, sigma = 0, 0.1  # mean and standard deviation
            pass

        self.erg = erg
        self._erg2speed()

    def move_ptcl(self, dt, EF):
        """Move each partile in a length of delta_L along its v-vector."""
        self.posn += self.uvec*self.speed*dt
        self.accl = EF*(self.charge*UNIT_CHARGE)/(self.mass*AMU)
        self.vel += self.accl*dt
        self.speed = sqrt(self.vel[0]**2 + self.vel[1]**2)
        self.uvec = self.vel/self.speed
        self._speed2erg()
        self._uvec2ang()
        self.step += 1

    def check_bndy(self, domain, imode=None):
        """
        Check the b.c. for moving ptcl.

        make the ptcl dead if it gets beyond the top bdry
        three modes for vertical bdry are available:
            lost, periodic and reflective
        left bdry is alway at 0.0
        """
        bottom, top, left, right = domain
        x, z = self.posn
        if imode == 'Lost':
            if not (left < x < right):
                self.isAlive = False
        elif imode == 'Periodic':
            self.posn[0] = left + ((x - left) % (right - left))
        elif imode == 'Reflective':
            pass
        else:
            return f'Error: imode={imode} is not found!'
        
    def check_wafer(self, wafer_loc):
        """Check if hit the wafer surface."""
        x, z = self.posn
        if z <= wafer_loc:
            self.isAlive = False
            self.cllct.append([self.erg, self.ang])