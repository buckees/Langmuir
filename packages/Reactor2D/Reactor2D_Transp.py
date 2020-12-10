"""
2D Plasma Transport Module.

Transp_2d contains:
    Diffusion only
    Ambipolar
    Drfit-Diffusion
    Momentum Solver

    Continuity Eq. dn/dt = -dF/dx * dt + S * det
    Input: depends on transport mode
    Output: dF/dx for continuity equation.
"""

import numpy as np
from copy import deepcopy

from packages.Constants import (UNIT_CHARGE, KB_EV)

class TRANSP2D(object):
    """Define the base tranport module/object."""
    
    def __init__(self, name='Txp2d'):
        """
        Init TRANSP2D.
        
        name: str, var, name of the TRANSP2D.
        """
        self.name = name
    
    def from_PLASMA(self, PLA):
        """Copy var from PLASMA2D."""
        self.ne, self.ni = deepcopy(PLA.ne), deepcopy(PLA.ni)
        self.pot = deepcopy(PLA.pot)
        self.Se, self.Si = deepcopy(PLA.Se), deepcopy(PLA.Si)
        self.Ex, self.Ez = deepcopy(PLA.Ex), deepcopy(PLA.Ez)
        self.fluxex = np.zeros_like(PLA.ne)
        self.fluxez = np.zeros_like(PLA.ne)
        self.fluxix = np.zeros_like(PLA.ne)
        self.fluxiz = np.zeros_like(PLA.ne)
        self.dfluxe = np.zeros_like(PLA.ne)
        self.dfluxi = np.zeros_like(PLA.ne)

        self._calc_txp_coeff(PLA)
    
    def to_PLASMA(self, PLA):
        """Copy var to PLASMA2D."""
        PLA.ne, PLA.ni = deepcopy(self.ne), deepcopy(self.ni)
        PLA.pot = deepcopy(self.pot)
        PLA.Ex, PLA.Ez = deepcopy(self.Ex), deepcopy(self.Ez)
        PLA.fluxex = deepcopy(self.fluxex)
        PLA.fluxez = deepcopy(self.fluxez)
        PLA.fluxix = deepcopy(self.fluxix)
        PLA.fluxiz = deepcopy(self.fluxiz)
        PLA.dfluxe = deepcopy(self.dfluxe)
        PLA.dfluxi = deepcopy(self.dfluxi)

    def _calc_txp_coeff(self, PLA):
        """
        Calc diffusion coefficient and mobility.

        PLA: PLASMA2D object/class
        De,i: m^2/s, (nz, nx) matrix, D = k*T/(m*coll_m)
        Mue,i: m^2/(V*s), (nz, nx) matrix, Mu = q/(m*coll_m)
        D and Mu depend only on PLA.
        """
        # calc diff coeff: D = k*T/(m*coll_m)
        self.De = np.divide(KB_EV*PLA.Te, PLA._Me*PLA.coll_em)  
        self.Di = np.divide(KB_EV*PLA.Ti, PLA._Mi*PLA.coll_im)  
        # calc mobility: Mu = q/(m*coll_m)
        self.Mue = UNIT_CHARGE/PLA._Me/PLA.coll_em
        self.Mui = UNIT_CHARGE/PLA._Mi/PLA.coll_im
    
    def solve_fluid(self, dt):
        """
        Evolve the density.
        
        SRC: REACT2D object/class
        dt: s, float, timestep
        """
        self.ne += (-self.dfluxe + self.Se)*dt
        self.ni += (-self.dfluxi + self.Si)*dt

class DIFF2D(TRANSP2D):
    """
    Calc the dflux for Diffusion Only Module.
    
    dn/dt = -D * d2n/dx2 + Se
    D: m^2/s, diffusion coefficient is calc from Tranps_1d
    Output: D * d2n/dx2
    """
    
    def calc_diff(self, PLA):
        """Calc diffusion term: D * d2n/dx2 and diffusion flux D * dn/dx."""
        # Calc transp coeff first
        self.calc_txp_coeff(PLA)
        # Calc flux
        self.fluxex, self.fluxez = -self.De * PLA.mesh.cnt_diff(self.ne)
        self.fluxix, self.fluxiz = -self.Di * PLA.mesh.cnt_diff(self.ni)
        # Calc dflux
        self.dfluxe = -self.De * PLA.mesh.cnt_diff_2nd(self.ne)
        self.dfluxi = -self.Di * PLA.mesh.cnt_diff_2nd(self.ni)

    
class AMBI2D(TRANSP2D):
    """
    Calc the dflux for Ambipolar Diffusion Module.

    dn/dt = -Da * d2n/dx2 + Se
    Di: m^2/s, ion diffusion coefficient is calc from Tranps_1d
    Da = Di(1 + Te/Ti).
    Output: Da * d2n/dx2 and E-field
    """

    def calc_ambi(self, PLA):
        """
        Calc ambipolar diffusion coefficient.

        The ambipolar diffusion assumptions:
            1. steady state, dne/dt = 0. it cannot be used to
            describe PLAsma decay.
            2. ni is calculated from continuity equation.
            3. PLAsma is charge neutral, ne = ni
            4. Ionization Se is needed to balance diffusion loss.
        Da = (De*Mui + Di*Mue)/(Mue + Mui)
        Da = Di(1 + Te/Ti).
        Ea = (Di - De)/(Mui + Mue)*dn/dx/n
        Orginal Ambipolar Coeff Da = (De*Mui + Di*Mue)/(Mue + Mui)
        self.Da = (PLAsma_1d.De*PLAsma_1d.Mui + PLAsma_1d.Di*PLAsma_1d.Mue) / \
                  (PLAsma_1d.Mue + PLAsma_1d.Mui)
        Assume Te >> Ti, Ambipolar Coeff can be simplified as
        Da = Di(1 + Te/Ti).
        """
        # Calc transp coeff first
        self._calc_txp_coeff(PLA)
        # Calc ambi coeff
        self.Da = self.Di*(1.0 + np.divide(PLA.Te, PLA.Ti))
        dnix, dniz = PLA.mesh.cnt_diff(self.ni)
        self.Ex = np.divide(self.Di - self.De, self.Mui + self.Mue)
        self.Ez = deepcopy(self.Ex)
        self.Ex *= np.divide(dnix, self.ni)
        self.Ez *= np.divide(dniz, self.ni)
        # # Calc flux
        self.fluxex, self.fluxez = -self.Da * PLA.mesh.cnt_diff(self.ne)
        self.fluxix, self.fluxiz = -self.Da * PLA.mesh.cnt_diff(self.ni)
        # Calc dflux
        self.dfluxe = -self.Da * PLA.mesh.cnt_diff_2nd(self.ne)
        self.dfluxi = -self.Da * PLA.mesh.cnt_diff_2nd(self.ni)
        
