"""
2D Plasma Electron Energy Module.

EERGY2D contains:
    Electron energy equation
    d(3/2nekTe)/dt = -dQ/dx + Power_in(ext.) - Power_loss(react)
    Input: ne, Te from Plasma1d, E_ext from field solver
    Output: Te
"""

import numpy as np
from copy import copy, deepcopy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
colMap = copy(cm.get_cmap("jet"))
colMap.set_under(color='white')

from packages.Constants import KB_EV

class EERGY2D(object):
    """Define the Eon Energy Module."""
    
    def __init__(self, name='Eergy2d'):
        """
        Init the EERGY2D.
        
        name: str, var, name of the EERGY2D.
        """
        self.name = name
    
    def from_PLASMA(self, PLA):
        """Copy var from PLASMA2D."""
        self.Te = deepcopy(PLA.Te)
        # eon energy = 3/2 * ne * kTe
        self.ergy_e = 1.5*KB_EV*np.multiply(PLA.ne, PLA.Te)
        self.fluxex = deepcopy(PLA.fluxex)
        self.fluxez = deepcopy(PLA.fluxez)
        self.dfluxe = deepcopy(PLA.dfluxe)
        
    def to_PLASMA(self, PLA):
        """Copy var to PLASMA2D."""
        PLA.Te = deepcopy(self.Te)
    
    def _calc_th_cond_coeff(self, PLA):
        """
        Calc thermal conduction coefficient.

        PLA: PLASMA2D object/class.
        th_cond_e: W/m/K, (nz, nx) matrix, heat conductivity for eon
        """
        # calc thermal conductivity for eon
        self.th_cond_e = np.ones_like(PLA.ne)*1e-3
        self._set_nonPlasma(PLA)

    def _set_nonPlasma(self, PLA):
        """Impose fixed th_cond_coeff on the non-PLAsma materials."""
        for idx, mat in np.ndenumerate(PLA.mesh.mat):
            if mat:
                self.th_cond_e[idx] = 1e-3
                self.Te[idx] = 0.1
    
    def _calc_th_flux(self, PLA):
        """
        Calc eon thermal flux, Qe.
        
        Qe = 5/2kTe * fluxe - ke * dTe/dx
        dQe = 5/2kTe * dfluxe - ke * d2Te/dx2
        PLA: PLASMA2D obj/class
        """
        # calc convection term
        self.Qex = 2.5*KB_EV*np.multiply(self.Te, self.fluxex)
        self.Qez = 2.5*KB_EV*np.multiply(self.Te, self.fluxez)
        self.dQe = 2.5*KB_EV*np.multiply(self.Te, self.dfluxe)
        # calc conduction term
        self.dTex, self.dTez = PLA.mesh.cnt_diff(self.Te)
        self.d2Te = PLA.mesh.cnt_diff_2nd(self.Te)
        self.Qex -= np.multiply(self.th_cond_e, self.dTex)
        self.Qez -= np.multiply(self.th_cond_e, self.dTez)
        self.dQe -= np.multiply(self.th_cond_e, self.d2Te)

    def _limit_Te(self, T_min=0.001, T_max=100.0):
        """Limit Te in the PLAsma."""
        self.Te = np.clip(self.Te, T_min, T_max)
        
    def solve_Te(self, PLA, dt):
        """
        Solve for Te.
        
        dt: s, var, time step for explict method
        PLA: PLASMA2D object/class.
        TXP: TRANSP2D object/class.
        """
        self._calc_th_cond_coeff(PLA)
        self._calc_th_flux(PLA)
        self.ergy_e += (-self.dQe + PLA.pwr_in)*dt
        self.Te = np.divide(self.ergy_e, PLA.ne)/1.5/KB_EV
        self._set_nonPlasma(PLA)
        self._limit_Te()

    def plot_dQe(self, PLA, figsize=(8, 8), ihoriz=1, 
                    dpi=300, fname='dQe.png', imode='Contour'):
        """
        Plot power vs. position.
            
        var include input power and total power.
        figsize: a.u., (2, ) tuple, size of fig
        ihoriz: a.u., var, 0 or 1, set the layout of fig horizontal or not
        dpi: a.u., dots per inch
        fname: str, var, name of png file to save
        imode: str, var, ['Contour', 'Scatter']
        """
        _x, _z = PLA.mesh.x, PLA.mesh.z
        if ihoriz:
            fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        else:
            fig, axes = plt.subplots(2, 1, figsize=figsize, dpi=dpi,
                                     constrained_layout=True)
        
        # plot densities
        if imode == 'Contour':
            for _ax, _den, _title in zip(axes, (self.dQe, self.pwr), 
                                ('Power Loss', 'Power to Eon')):
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
