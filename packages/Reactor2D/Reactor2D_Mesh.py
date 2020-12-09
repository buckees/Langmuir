"""
Mesh Module.

Create standalone mesh or
Create mesh for input geometry.
"""

import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

class MESH2D():
    """Define 2d Mesh."""

    def __init__(self, name='Mesh2d'):
        """
        Init the Shape.
        
        name: str, var, name of the Mesh2d.
        """
        self.name = name

    def import_geom(self, geom):
            """Import geometry."""
            self.geom = geom

    def generate_mesh(self, ngrid=(11, 11)):
        """Generate mesh according to the imported geometry."""
        self.width, self.height = self.geom.domain
        self.ngrid = np.asarray(ngrid)
        self.nx, self.nz = self.ngrid
        self.res = np.divide(self.geom.domain, self.ngrid - 1)
        self.dx, self.dz = self.res
        tempx = np.linspace(self.geom.bl[0], self.geom.bl[0] + self.width, 
                            self.nx)
        tempz = np.linspace(self.geom.bl[1], self.geom.bl[1] + self.height, 
                            self.nz)
        self.x, self.z = np.meshgrid(tempx, tempz)
        self.mat = np.zeros_like(self.x)
        self._find_bndy()
        self._assign_mat()
        self._calc_plasma_area()
       
    def create_mesh(self, bl=(0.0, 0.0), domain=(1.0, 1.0), ngrid=(11, 11)):
        """Create standalone mesh."""
        self.bl = np.asarray(bl)
        self.domain = np.asarray(domain)
        self.ngrid = np.asarray(ngrid)
        self.res = np.divide(self.domain, self.ngrid - 1)
        self.width, self.height = self.domain
        self.nx, self.nz = self.ngrid
        self.delx, self.delz = self.res
        tempx = np.linspace(0.0, self.width, self.nx)
        tempz = np.linspace(0.0, self.height, self.nz)
        self.x, self.z = np.meshgrid(tempx, tempz)
        self._find_bndy()

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

    def _assign_mat(self):
        """Assign materials to nodes."""
        for _idx, _x in np.ndenumerate(self.x):
            _z = self.z[_idx]
            _posn = np.array([_x, _z])
            _label, self.mat[_idx] = self.geom.get_label(_posn)
    
    def _calc_plasma_area(self):
        """Calc the total area of plasma region."""
        self.area = 0
        for _idx, _mat in np.ndenumerate(self.mat):
            if not _mat:
                self.area += self.dx * self.dz

    def plot(self, figsize=(8, 8), dpi=600, fname='Mesh.png',  ihoriz=1):
        """Plot mesh."""
        colMap = plt.get_cmap('Set1')
        
        if ihoriz:
            fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        ax = axes[0]
        ax.scatter(self.x, self.z, c=self.mat, s=10, cmap=colMap)
        ax = axes[1]
        ax.scatter(self.x, self.z, c=self.bndy, s=10, cmap=colMap)
        fig.savefig(fname, dpi=dpi)
        plt.close()

    def cnt_diff(self, f):
        """
        Caculate dy/dx using central differencing.
        
        input: y
        dy/dx = (y[i+1] - y[i-1])/(2.0*dx)
        dy[0] = dy[1]; dy[-1] = dy[-2]
        output: dy
        """
        dfx = np.zeros_like(self.x)
        dfz = np.zeros_like(self.z)
        # Although dy[0] and dy[-1] are signed here,
        # they are eventually specified in boundary conditions
        # dy[0] = dy[1]; dy[-1] = dy[-2]
        for i in range(1, self.nx-1):
            dfx[:, i] = (f[:, i+1] - f[:, i-1])/self.dx/2.0
        for j in range(1, self.nz-1):
            dfz[j, :] = (f[j+1, :] - f[j-1, :])/self.dz/2.0
        dfx[:, 0], dfx[:, -1] = deepcopy(dfx[:, 1]), deepcopy(dfx[:, -2])
        dfz[0, :], dfz[-1, :] = deepcopy(dfz[1, :]), deepcopy(dfz[-2, :])
        return dfx, dfz
    
    def cnt_diff_2nd(self, f):
        """
        Caculate d2y/dx2 using 2nd order central differencing.

        input: y
        d2y/dx2 = (y[i+1] - 2 * y[i] + y[i-1])/dx^2
        d2y[0] = d2y[1]; d2y[-1] = d2y[-2]
        output: d2y/dx2
        """
        d2fx = np.zeros_like(self.x)
        d2fz = np.zeros_like(self.z)
        # Although dy[0] and dy[-1] are signed here,
        # they are eventually specified in boundary conditions
        # d2y[0] = d2y[1]; d2y[-1] = d2y[-2]
        for i in range(1, self.nx-1):
            d2fx[:, i] = (f[:, i+1] - 2 * f[:, i] + f[:, i-1])/self.dx**2
        for j in range(1, self.nz-1):
            d2fz[j, :] = (f[j+1, :] - 2 * f[j, :] + f[j-1, :])/self.dz**2
        d2fx[:, 0], d2fx[:, -1] = deepcopy(d2fx[:, 1]), deepcopy(d2fx[:, -2])
        d2fz[0, :], d2fz[-1, :] = deepcopy(d2fz[1, :]), deepcopy(d2fz[-2, :])
        d2f = d2fx + d2fz
        return d2f
