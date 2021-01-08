"""
Mesh Module.

Create standalone mesh or
Create mesh for input geometry.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from copy import copy
colMap = copy(cm.get_cmap("jet"))
colMap.set_under(color='white')

from packages.Mesh.Mesh import MESH2D

class MESH2D(MESH2D):
    """Define 2d Mesh."""

    def __init__(self, name='Mesh2d'):
        """
        Init the Shape.
        
        name: str, var, name of the Mesh2d.
        """
        self.name = name

    def readin_mesh(self, fname='restart'):
        """Read in mat, x and z from npy files."""
        temp = np.load(fname+'.npz', allow_pickle=True)
        self.x = temp['x']
        self.z = temp['z']
        self.mat = temp['mat']
        self.res = temp['res']
        self.ngrid = temp['ngrid']
        self.bl = temp['bl']
        self.tr = temp['tr']
        self.mat_dict = temp['mat_dict'].item()
        self.nx, self.nz = self.ngrid
        self.dx, self.dz = self.res
        self.left, self.bottom = self.bl
        self.right, self.top = self.tr
        self.width = self.right - self.left
        self.height = self.top - self.bottom
        # find the surf nodes and surf_vac nodes
        self._find_bndy()
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

    def _calc_plasma_area(self):
        """Calc the total area of plasma region."""
        self.isPlasma = (self.mat == 0).astype(int)
        self.area = self.dx * self.dz * self.isPlasma.sum()
        
    def plot_var(self, var, var_name,
                 fname='Plasma.png',figsize=(16, 8), ihoriz=1, dpi=300, 
                 imode='Contour', iplot_geom=0):
        """
        Plot plasma variables vs. position.
            
        var: list of var, such as [ne, ni]
        var_name: list of str, such as ['E Density', 'Ion Density']
        fname: str, var, name of png file to save
        figsize: a.u., (2, ) tuple, size of fig
        ihoriz: a.u., var, 0 or 1, set the layout of fig horizontal or not
        dpi: a.u., dots per inch
        imode: str, var, ['Contour', 'Scatter']
        iplot_geom: int, var, control whether to plot geom
        """
        nvar = len(var)
        if ihoriz:
            fig, axes = plt.subplots(1, nvar, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(nvar, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        # plot var
        for ax, den, title in zip(axes, var, var_name):
            if imode == 'Contour':
                cs = ax.contourf(self.x, self.z, den, cmap=colMap)
            elif imode == 'Scatter':
                cs = ax.scatter(self.x, self.z, c=den, cmap=colMap)
            ax.set_title(title)
            fig.colorbar(cs, ax=ax, shrink=0.9)
            ax.set_xlabel('Position (m)')
            ax.set_ylabel('Height (m)')
            ax.set_aspect('equal')
        # save and close()
        fig.savefig(fname, dpi=dpi)
        plt.close()
