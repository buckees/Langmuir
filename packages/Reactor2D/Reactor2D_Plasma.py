"""
2D Plasma Module - Data Transfer/Communication Hub.

PLASMA2D contains:
    mesh
    all updated variables
    plot
    diagnostics
"""

import numpy as np

from packages.Constants import (AMU, UNIT_CHARGE, EON_MASS)

class PLASMA2D(object):
    """Define PLASMA2D."""

    def __init__(self, name='Plasma2d'):
        """
        Init the PLASMA2D.
        
        name: str, var, name of the PLASMA2D.
        """
        self.name = name

    def import_mesh(self, mesh):
        """Import mesh."""
        self.mesh = mesh

    def init_plasma(self, ne=1e17, press=10, 
                    Te=1, Ti=0.1, Mi=40, wrf=13.56e6):
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
        self._wrf = wrf # rf frequency in angular frequency,  rad/s
        # public attributes
        x = self.mesh.x
        self.ne = np.ones_like(x)*ne  # init uniform ne on 1d mesh
        self.ni = np.ones_like(x)*ne  # init ni to neutralize ne
        self.nn = np.ones_like(x)*(press*3.3e19)  # init neutral density
        self.Te = np.ones_like(x)*Te  # init eon temperature
        self.Ti = np.ones_like(x)*Ti  # init ion temperature
        self.coll_em = np.ones_like(x)*(press/10.0*1e7)  # eon coll freq (mom)
        self.coll_im = np.ones_like(x)*(press/10.0*1e7)  # ion coll freq (mom)
        self.pot = np.zeros_like(x)  # initial uniform potential
        self.Ex = np.zeros_like(x)  # initial uniform E-field
        self.Ez = np.zeros_like(x)  # initial uniform E-field
        self.Ey = np.zeros_like(x)  # initial uniform E-field
        self.conde = np.zeros_like(x)
        self.pwr_in = np.zeros_like(x)
        self.fluxex, self.fluxez = np.zeros_like(x), np.zeros_like(x)
        self.fluxix, self.fluxiz = np.zeros_like(x), np.zeros_like(x)
        self.dfluxe, self.dfluxi = np.zeros_like(x), np.zeros_like(x)
        # modify init
        self.update_plasma()

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
    
    def update_plasma(self):
        """
        Make sure the all attr in PLASMA2D are updated.
        
        aaa
        """
        self._calc_conde()
        self._calc_pwr_in()
        self._set_bc()
        self._set_nonPlasma()
        self._limit_plasma()

    def _calc_pwr_in(self):
        """Calc input power due to E-field."""
        EF2 = self.Ex**2 + self.Ez**2 + self.Ey**2
        self.pwr_in = np.multiply(self.conde, EF2)
    
    def get_eps(self):
        """Get epsilon from materials."""
        self.eps = np.ones_like(self.ne)
        for idx, mat in np.ndenumerate(self.mesh.mat):
            if mat == 2:
                self.eps[idx] = 3.8
        
    def _calc_conde(self):
        """Calc eon conductivity."""
        temp_conde = UNIT_CHARGE**2/EON_MASS
        temp_conde *= self.ne
        self.conde = np.zeros_like(self.ne, dtype=complex)
        self.conde += self.coll_em
        self.conde += complex(0.0, self._wrf)
        self.conde = np.divide(temp_conde, self.conde)