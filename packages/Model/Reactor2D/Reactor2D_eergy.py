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
        self.dQe = 0.0
        # self.dQe = PLA.mesh.cnt_diff((self.Qex, self.Qez), imode='Vector')
        # self.dQe = 2.5*KB_EV*np.multiply(self.Te, self.dfluxe)
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