"""
Mesh Module.

Create mesh for a given geometry.
"""

import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

class Mesh():
    """Define all shared basic properties."""

    def __init__(self, import_geom):
        """Import geometry as init."""
        self.geom = import_geom
        self.name = import_geom.name + '_Mesh'
            
    def __str__(self):
        """Print out info."""
        res = f'This is a {self.geom.dim}D mesh called {self.geom.name},'
        return res
    
        
class Mesh2D(Mesh):
    """Define 2D mesh."""

    def gen_mesh(self, ngrid=(10, 10)):
        """Generate mesh according to the imported geometry."""
        self.ngrid = np.asarray(ngrid)
        self.nx, self.nz = self.ngrid
        self.res = np.divide(self.geom.domain.domain, self.ngrid)
        self.dx, self.dz = self.res
        tempx = np.linspace(self.geom.domain.bl[0] + 0.5*self.dx, 
                            self.geom.domain.tr[0] - 0.5*self.dx, 
                            self.nx)
        tempz = np.linspace(self.geom.domain.bl[1] + 0.5*self.dz, 
                            self.geom.domain.tr[1] - 0.5*self.dz, 
                            self.nz)
        self.x, self.z = np.meshgrid(tempx, tempz)
        self.mat = np.zeros_like(self.x)
        self._find_bndy()
        self._assign_mat()
        self._calc_plasma_area()
       
    def _find_bndy(self):
        """Add boundaries."""
        self.bndy = np.zeros_like(self.x)
        self.bndy_list = list()
        for i in range(self.nx-1):
            self.bndy_list.append((0, i))
        for j in range(self.nz-1):
            self.bndy_list.append((j, self.nx-1))
        for i in reversed(range(1, self.nx)):
            self.bndy_list.append((self.nz-1, i))
        for j in reversed(range(1, self.nz)):
            self.bndy_list.append((j, 0))
        # sign value at bndy as 1
        for idx in self.bndy_list:
            self.bndy[idx] = 1

    def _assign_mat(self):
        """Assign materials to nodes."""
        for idx, x in np.ndenumerate(self.x):
            z = self.z[idx]
            posn = np.array([x, z])
            mater = self.geom.get_mater(posn)
            self.mat[idx] = self.geom.mat_dict[mater]
    
    def _calc_plasma_area(self):
        """Calc the total area of plasma region."""
        self.area = 0
        for idx, mat in np.ndenumerate(self.mat):
            if not mat:
                self.area += self.dx * self.dz

    def plot(self, figsize=(8, 8), dpi=600, ihoriz=1, s_size=10):
        """Plot mesh."""
        colMap = plt.get_cmap('Set1')
        
        if ihoriz:
            fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        ax = axes[0]
        ax.scatter(self.x, self.z, c=self.mat, s=s_size, cmap=colMap)
        ax = axes[1]
        ax.scatter(self.x, self.z, c=self.bndy, s=s_size, cmap=colMap)
        fig.savefig(self.name, dpi=dpi)
        plt.close()
        

class Mesh1D(Mesh):
    """Define 1D Mesh."""

    def gen_mesh(self, nx=10):
        """Generate mesh according to the imported geometry."""
        self.nx = nx
        self.x, self.dx = np.linspace(self.geom.domain.domain[0], 
                                      self.geom.domain.domain[1], 
                                      self.nx, retstep=True)
        self.mat = np.zeros_like(self.x)
        self._find_bndy()
        self._assign_mat()
        self._calc_plasma_area()

    def _find_bndy(self):
        """Add boundaries."""
        self.bndy = np.zeros_like(self.x)
        self.bndy[0], self.bndy[-1] = 1, 1        

    def _assign_mat(self):
        """Assign materials to nodes."""
        for idx, x in np.ndenumerate(self.x):
            mater = self.geom.get_mater(x)
            self.mat[idx] = self.geom.mater_dict[mater]
    
    def _calc_plasma_area(self):
        """Calc the total area of plasma region."""
        self.area = 0
        for mat in self.mat:
            if not mat:
                self.area += self.dx

    def plot(self, figsize=(8, 8), dpi=600, ihoriz=1):
        """Plot mesh."""
        colMap = plt.get_cmap('Set1')
        
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi,
                                 constrained_layout=True)

        ax.scatter(self.x, np.zeros_like(self.x), 
                   c=self.mat, s=10, cmap=colMap)
        fig.savefig(self.name, dpi=dpi)
        plt.close()
