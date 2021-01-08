"""
Mesh Module.

Create standalone mesh or
Create mesh for input geometry.
"""

import numpy as np

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
