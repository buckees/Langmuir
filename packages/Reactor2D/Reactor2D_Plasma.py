"""
2D Plasma Module - Data Transfer/Communication Hub.

PLASMA2D contains:
    mesh
    all updated variables
    plot
    diagnostics
"""

import numpy as np
from copy import copy, deepcopy
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import matplotlib.cm as cm
colMap = copy(cm.get_cmap("jet"))
colMap.set_under(color='white')

from packages.Constants import (AMU, UNIT_CHARGE, EON_MASS, COLOR_DICT)

class PLASMA2D(object):
    """Define PLASMA2D."""

    def __init__(self, name='Plasma2d'):
        """
        Init the Shape.
        
        name: str, var, name of the Mesh2d.
        """
        self.name = name

    def import_mesh(self, mesh):
        """Import mesh."""
        self.mesh = mesh

    def init_plasma(self, ne=1e17, press=10, Te=1, Ti=0.1, Mi=40):
        """
        Initiate plasma attributes.

        ne: 1/m^3, eon denisty
        ni: 1/m^3, ion density = eon density initially
        press: mTorr, pressure
        nn: 1/m^3, neutral density determined by pressure
            At 1 atm, number density = 0.025e27 m^-3.
            At 1 Torr, number density = 3.3e22 m^-3.
            At 1 mTorr, number density = 3.3e19 m^-3.
        Te: eV, eon temperature
        Ti: eV, ion temperature
        coll_e,im: 1/s, coll freq (momentum)
                coll_e,im = 1e7 at 10 mTorr
                            1e8 at 100 mTorr
                            1e9 at 1000 mTorr
        Mi: kg, ion mass
        """
        # private attributes
        self._press = press # pressure in Torr
        self._Me = EON_MASS # eon mass in kg
        self._Mi = Mi*AMU # ion mass in kg
        # public attributes
        x = self.mesh.x
        self.ne = np.ones_like(x)*ne  # init uniform ne on 1d mesh
        self.ni = np.ones_like(x)*ne  # init ni to neutralize ne
        self.nn = np.ones_like(x)*(press*3.3e19)  # init neutral density
        self.Te = np.ones_like(x)*Te  # init eon temperature
        self.Ti = np.ones_like(x)*Ti  # init ion temperature
        self.coll_em = np.ones_like(x)*(press/10.0*1e7)  # eon coll freq (mom)
        self.coll_im = np.ones_like(x)*(press/10.0*1e7)  # ion coll freq (mom)
        
        self._set_bc()
        self._set_nonPlasma()
        self._limit_plasma()

    def _set_bc(self):
        """Impose b.c. on the plasma."""
        for idx in self.mesh.bndy_list:
            self.ne[idx] = 1e11
            self.ni[idx] = 1e11
            self.nn[idx] = 1e11
            self.Te[idx] = 0.1
            self.Ti[idx] = 0.01

    def _set_nonPlasma(self):
        """Impose fixed values on the non-plasma materials."""
        for idx, mat in np.ndenumerate(self.mesh.mat):
            if mat:
                self.ne[idx] = 1e11
                self.ni[idx] = 1e11
                self.nn[idx] = 1e11
                self.Te[idx] = 0.1
                self.Ti[idx] = 0.01

    def _limit_plasma(self, n_min=1e11, n_max=1e22, T_min=0.001, T_max=100.0):
        """Limit variables in the plasma."""
        self.ne = np.clip(self.ne, n_min, n_max)
        self.ni = np.clip(self.ni, n_min, n_max)
        self.nn = np.clip(self.nn, n_min, n_max)
        self.Te = np.clip(self.Te, T_min, T_max)
        self.Ti = np.clip(self.Ti, T_min, T_max)

    def plot_plasma(self, figsize=(8, 8), ihoriz=1, 
                    dpi=300, fname='Plasma.png', imode='Contour',
                    iplot_geom=0):
        """
        Plot plasma variables vs. position.
            
        var include density, temperature.
        figsize: a.u., (2, ) tuple, size of fig
        ihoriz: a.u., var, 0 or 1, set the layout of fig horizontal or not
        dpi: a.u., dots per inch
        fname: str, var, name of png file to save
        imode: str, var, ['Contour', 'Scatter']
        """
        _x, _z = self.mesh.x, self.mesh.z
        if ihoriz:
            fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        
        # plot densities
        if imode == 'Contour':
            for _ax, _den, _title in zip(axes, (self.ne, self.ni), 
                                         ('E Density', 'Ion Density')):
                _cs = _ax.contourf(_x, _z, _den, cmap=colMap, vmin=1.1e11)
                _ax.set_title(_title)
                fig.colorbar(_cs, ax=_ax, shrink=0.9)
            
        elif imode == 'Scatter':
            for _ax, _den, _title in zip(axes, (self.ne, self.ni), 
                                         ('E Density', 'Ion Density')):
                _ax.scatter(_x, _z, c=_den, s=5, cmap=colMap, vmin=1.1e11)
                _ax.set_title(_title)
            
        for ax in axes:
            ax.set_xlabel('Position (m)')
            ax.set_ylabel('Height (m)')
            ax.set_aspect('equal')
            
            # add geom plot
            if iplot_geom:
                for shape in self.mesh.geom.sequence:
                    if shape.type == 'Rectangle':
                        temp_col = COLOR_DICT[self.mesh.geom.label[shape.label]]
                        ax.add_patch(
                            patch.Rectangle(shape.bl, shape.width, shape.height,
                                            facecolor=temp_col, edgecolor='w')
                            )
        
        fig.savefig(fname, dpi=dpi)
        plt.close()
    
    def get_Te(self, een):
        """
        Get Te from Eergy2d().
        
        een: Eergy2d() boject.
        """
        self.Te = deepcopy(een.Te)
    
    def plot_Te(self, figsize=(8, 8), ihoriz=1, 
                    dpi=300, fname='Te.png', imode='Contour'):
        """
        Plot plasma variables vs. position.
            
        var include density, temperature.
        figsize: a.u., (2, ) tuple, size of fig
        ihoriz: a.u., var, 0 or 1, set the layout of fig horizontal or not
        dpi: a.u., dots per inch
        fname: str, var, name of png file to save
        imode: str, var, ['Contour', 'Scatter']
        """
        _x, _z = self.mesh.x, self.mesh.z
        if ihoriz:
            fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        
        # plot densities
        if imode == 'Contour':
            for _ax, _den, _title in zip(axes, (self.Te, self.Ti), 
                                ('E Temperature', 'Ion Temperature')):
                _cs = _ax.contourf(_x, _z, _den, cmap=colMap)
                _ax.set_title(_title)
                fig.colorbar(_cs, ax=_ax, shrink=0.9)
            
        elif imode == 'Scatter':
            for _ax, _den, _title in zip(axes, (self.Te, self.Ti), 
                                ('E Temperature', 'Ion Temperature')):
                _ax.scatter(_x, _z, c=_den, s=5, cmap=colMap)
                _ax.set_title(_title)
            
        for ax in axes:
            ax.set_xlabel('Position (m)')
            ax.set_ylabel('Height (m)')
            ax.set_aspect('equal')
        fig.savefig(fname, dpi=dpi)
        plt.close()

    def init_pot(self, phi=0.0):
        """Initiate potential attributes."""
        nx = self.mesh.nx
        self.pot = np.ones(nx)*phi  # initial uniform potential
        self.ef = np.zeros_like(self.pot)  # initial uniform E-field
        self.ef_ambi = np.zeros_like(self.pot)  # initial ambipolar E-field

    def plot_pot(self):
        """Plot potential, E-field."""
        x = self.mesh.x
        fig, axes = plt.subplots(1, 2, figsize=(8, 4),
                                 constrained_layout=True)
        # plot potential
        ax = axes[0]
        ax.plot(x, self.pot, 'y-')
        ax.legend(['Potential'])
        # plot E-field
        ax = axes[1]
        ax.plot(x, self.ef, 'g-')
        ax.legend(['E-field'])
        plt.show()

    def den_evolve(self, delt, txp, src):
        """
        Evolve the density in Plasma by solving the continuity equation.

        dn/dt = -dFlux/dx + Se
        dn(t + dt) = dn(t) - dFlux/dx*dt + Se*dt
        delt: s, var, time step for explict method
        txp: Transp2d() object
        src: React2d() object
        """
        self.ne += (-txp.dfluxe + src.Se)*delt
        self.ni += (-txp.dfluxi + src.Si)*delt
        self._set_bc()
        self._set_nonPlasma()
        self._limit_plasma()

    def readin_EF(self,fname):
        """Read in E-Field from external file."""
        temp_EF = np.fromfile(fname, dtype='f4')
        temp_EF = np.reshape(temp_EF, (41, 51))
        temp_EF = np.flip(temp_EF, 0)
        self.EF = temp_EF
    
    def get_eps(self):
        """Get epsilon from materials."""
        self.eps = np.ones_like(self.ne)
        for idx, mat in np.ndenumerate(self.mesh.mat):
            if mat == 2:
                self.eps[idx] = 3.8
        
    def plot_eps(self, figsize=(8, 8), ihoriz=1, 
                    dpi=300, fname='eps.png', imode='Contour'):
        """
        Plot eon conductivity vs. position.
            
        var include density, temperature.
        figsize: a.u., (2, ) tuple, size of fig
        ihoriz: a.u., var, 0 or 1, set the layout of fig horizontal or not
        dpi: a.u., dots per inch
        fname: str, var, name of png file to save
        imode: str, var, ['Contour', 'Scatter']
        """
        x, z = self.mesh.x, self.mesh.z
        if ihoriz:
            fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        
        # plot eon conductivity
        for ax, var, title in zip(axes, (self.eps, self.eps), 
                        ('Dielectric Constant', 'Dielectric Constant')):
            if imode == 'Contour':
                cs = ax.contourf(x, z, var, cmap=colMap)
            elif imode == 'Scatter':
                cs = ax.scatter(x, z, c=var, s=5, cmap=colMap)
            ax.set_title(title)
            fig.colorbar(cs, ax=ax, shrink=0.9)
        
        for ax in axes:
            ax.set_xlabel('Position (m)')
            ax.set_ylabel('Height (m)')
            ax.set_aspect('equal')
        fig.savefig(fname, dpi=dpi)
        plt.close()

    def calc_conde(self, w_rf):
        """Calc eon conductivity."""
        temp_conde = UNIT_CHARGE**2/EON_MASS
        temp_conde *= self.ne
        self.conde = np.zeros_like(self.ne, dtype=complex)
        self.conde += self.coll_em
        self.conde += complex(0.0, w_rf)
        self.conde = np.divide(temp_conde, self.conde)
    
    def plot_conde(self, figsize=(8, 8), ihoriz=1, 
                    dpi=300, fname='conde.png', imode='Contour'):
        """
        Plot eon conductivity vs. position.
            
        var include density, temperature.
        figsize: a.u., (2, ) tuple, size of fig
        ihoriz: a.u., var, 0 or 1, set the layout of fig horizontal or not
        dpi: a.u., dots per inch
        fname: str, var, name of png file to save
        imode: str, var, ['Contour', 'Scatter']
        """
        x, z = self.mesh.x, self.mesh.z
        if ihoriz:
            fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        
        # plot eon conductivity
        for ax, var, title in zip(axes, (self.conde.real, self.conde.imag), 
                        ('E conductivity REAL', 'E conductivity IMAG')):
            if imode == 'Contour':
                cs = ax.contourf(x, z, var, cmap=colMap)
            elif imode == 'Scatter':
                cs = ax.scatter(x, z, c=var, s=5, cmap=colMap)
            ax.set_title(title)
            fig.colorbar(cs, ax=ax, shrink=0.9)
        
        for ax in axes:
            ax.set_xlabel('Position (m)')
            ax.set_ylabel('Height (m)')
            ax.set_aspect('equal')
        fig.savefig(fname, dpi=dpi)
        plt.close()
