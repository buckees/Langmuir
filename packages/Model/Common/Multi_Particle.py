"""
Multi_Particle.py serves as a data center/hub,
    is shered by all particle-based models,
    uses DataFrame from pandas,
    supports an array of particles.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from packages.Constants import (PI, AMU, UNIT_CHARGE, EV2J, J2EV)

class MULTI_PARTICLE(object):
    """Create MULTI_PARTICLE() object."""

    def __init__(self, name='Multi_Particle'):
        """
        Init the MULTI_PARTICLE().
        
        name: str, name of the MULTI_PARTICLE.
        num: int, num of particles
        pname: list of str, name of particle
        ptype: list of str, type of particle, ['Eon', 'Ion', 'Neut']
        mass: arr of float, mass of particle, unit in AMU
        charge: arr of float, charge of particle
        isAlive: arr of bool, state of particle
        posn: arr of float, position of particle, unit in m, shape=(num, 3)
        vel: arr of float, velocity of particle, unit in m/s, shape=(num, 3)
        erg: arr of float, energy of particle, unit in eV
        """
        self.name = name
        self.num = 1
        self.pname = list()
        self.ptype = list()
        self.mass = np.array(list())
        self.charge = np.array(list())
        self.isAlive = np.array(list())
        self.posn = np.array(list())
        self.vel = np.array(list())
        self.erg = np.array(list())
        
    def add_PARTICLE(self, PARTICLE):
        """Add PARTICLE() object to MULTI_PARTICLE()."""
        pass
    
    def to_array(self, weight):
        """Convert data to Numpy array type."""
        pass
    
    def gen_particles(self, num, prop, posn, vel):
        """
        Generate multi particles of a single species directly.
        
        num: int, num of particles.
        prop: dict, properties of the particle.
        posn: numpy array, (num, 3), position of all particles.
        vel: numpy array, (num, 3), velocity of all particles
        """
        self.num = num
        self.pname = [prop['name'] for i in range(num)]
        self.ptype = [prop['type'] for i in range(num)]
        self.mass = np.ones(num)*prop['mass']
        self.charge = np.ones(num)*prop['charge']
        self.isAlive = np.ones(num, dtype=bool)
        self.posn = posn
        self.vel = vel
        return print(f'{num} of particles are generated!')
    
    def update_posn(self, posn):
        """Update positions for all particles."""
        self.posn = posn
    
    def update_vel(self, vel):
        """Update velocities for all particles."""
        self.vel = vel
    
    def update_state(self, state):
        """Update states for all particles."""
        self.isAlive = state
    
    def plot(self):
        """Plot basic statics."""
        fig, axes = plt.subplots(2, 3, figsize=(12, 8), dpi=600,
                       constrained_layout=True)
        variables = self.posn.T.tolist()
        variables += self.vel.T.tolist()
        titles = ['Position in X', 'Position in Z', 'Position in Y',
                  'Velocity in X', 'Velocity in Z', 'Velocity in Y']
        xlabels = ['Position (m)', 'Position (m)', 'Position (m)',
                  'Velocity (m/s)', 'Velocity (m/s)', 'Velocity (m/s)']
        for ax, var, title, xlabel in zip(axes.flatten(), 
                                          variables, titles, xlabels):
            ax.hist(var, density=True)
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel('Probability')
        fig.savefig(f'{self.name}.png', dpi=600)
        plt.close()

    def vel2erg(self):
        """Convert velocity to energy."""
        temp = np.power(self.vel, 2)
        temp = np.sum(temp, axis=1)
        self.erg = 0.5*(self.mass*AMU)*temp*J2EV

    def plot_erg(self):
        """Plot basic statics in energy."""
        self.vel2erg()
        fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=600,
                       constrained_layout=True)
        ax.hist(self.erg, density=True)
        ax.set_title('Energy Distribution')
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Probability')
        fig.savefig(f'{self.name}_erg.png', dpi=600)
        plt.close()
    
    def to_DataFrame(self):
        """Convert data to DataFrame type."""
        pass
