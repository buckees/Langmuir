"""
2D Plasma Module - Data Transfer/Communication Hub.

PLASMA2D contains:
    MESH
    all updated variables
    plot
    diagnostics
"""

import numpy as np

from packages.Constants import (PI, AMU, UNIT_CHARGE, EON_MASS, KB_EV)

class PLASMA2D(object):
    """Define PLASMA2D."""

    def __init__(self, name='Plasma2d'):
        """
        Init the PLASMA2D.
        
        name: str, var, name of the PLASMA2D.
        """
        self.name = name

    def init_plasma(self, MESH, ne=1e17, press=10, 
                    Te=1, Ti=0.1, Tn=0.025, Mi=40, freq=13.56e6):
        """
        Initiate plasma attributes.

        ne: 1/m^3, eon denisty
        ni: 1/m^3, ion density = eon density initially
        press: mTorr, pressure
        freq: Hz, frequency
        nn: 1/m^3, neutral density determined by pressure
            At 1 atm, number density = 0.025e27 m^-3.
            At 1 Torr, number density = 3.3e22 m^-3.
            At 1 mTorr, number density = 3.3e19 m^-3.
        Te, Ti, Tn: eV, eon, ion, neut temperature
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
        self._wrf = 2.0*PI*freq # rf frequency in angular frequency,  rad/s
        # essential attributes for restart
        x = MESH.x
        self.ne = np.ones_like(x)*ne  # init uniform ne on 1d MESH
        self.ni = np.ones_like(x)*ne  # init ni to neutralize ne
        self.nn = np.ones_like(x)*(press*3.3e19)  # init neutral density
        self.Te = np.ones_like(x)*Te  # init eon temperature
        self.Ti = np.ones_like(x)*Ti  # init ion temperature
        self.Tn = np.ones_like(x)*Tn  # init neut temperature
        self.pot = np.zeros_like(x)  # initial uniform potential
        self._init_nonessential(MESH)
        # modify init
        self.update_plasma(MESH)
    
    def _init_nonessential(self, MESH):
        x = MESH.x
        self.Se, self.Si = np.zeros_like(x), np.zeros_like(x)
        self.Ex, self.Ez = np.zeros_like(x), np.zeros_like(x)
        self.fluxex, self.fluxez = np.zeros_like(x), np.zeros_like(x)
        self.fluxix, self.fluxiz = np.zeros_like(x), np.zeros_like(x)
        self.dfluxe, self.dfluxi = np.zeros_like(x), np.zeros_like(x)
        
    def update_plasma(self, MESH):
        """
        Make sure the all attr in PLASMA2D are updated.
        
        aaa
        """
        self._set_nonPlasma(MESH)
        self._limit_plasma()
        self._calc_conde()
        self._calc_pwr_in()
        self._calc_coll()
        self._calc_txp_coeff()
        self._calc_ave(MESH)
        
    def _set_nonPlasma(self, MESH):
        """Impose fixed values on the non-plasma materials."""
        for idx, mat in np.ndenumerate(MESH.mat):
            if mat:
                self.ne[idx], self.ni[idx], self.nn[idx] = 1e11, 1e11, 1e11
                self.Te[idx], self.Ti[idx] = 0.1, 0.01
                self.Ex[idx], self.Ez[idx] = 0.0, 0.0

    def _limit_plasma(self, n_min=1e11, n_max=1e22, T_min=0.001, T_max=100.0):
        """Limit variables in the plasma."""
        self.ne = np.clip(self.ne, n_min, n_max)
        self.ni = np.clip(self.ni, n_min, n_max)
        self.nn = np.clip(self.nn, n_min, n_max)
        self.Te = np.clip(self.Te, T_min, T_max)
        self.Ti = np.clip(self.Ti, T_min, T_max)
        
    def _calc_conde(self):
        """Calc eon conductivity."""
        temp_conde = UNIT_CHARGE**2/EON_MASS
        temp_conde *= self.ne
        self.conde = np.zeros_like(self.ne, dtype=complex)
        self.conde += self.coll_em
        self.conde += complex(0.0, self._wrf)
        self.conde = np.divide(temp_conde, self.conde)
        
    def _calc_pwr_in(self):
        """Calc input power due to E-field."""
        EF2 = np.abs(self.Ex)**2 + np.abs(self.Ez)**2 + np.abs(self.Ey)**2
        self.pwr_in = np.multiply(np.real(self.conde), EF2)
    
    def _calc_txp_coeff(self):
        """
        Calc diffusion coefficient and mobility.

        MESH: MESH2D object/class
        De,i: m^2/s, (nz, nx) matrix, D = k*T/(m*coll_m)
        Mue,i: m^2/(V*s), (nz, nx) matrix, Mu = q/(m*coll_m)
        D and Mu depend only on PLA.
        """
        # calc diff coeff: D = k*T/(m*coll_m)
        self.De = np.divide(KB_EV*self.Te, self._Me*self.coll_em)  
        self.Di = np.divide(KB_EV*self.Ti, self._Mi*self.coll_im)  
        # calc mobility: Mu = q/(m*coll_m)
        self.Mue = UNIT_CHARGE/self._Me/self.coll_em
        self.Mui = UNIT_CHARGE/self._Mi/self.coll_im
    
    def _calc_ave(self, MESH):
        """Calc averaged variables."""
        sum_isPlasma = MESH.isPlasma.sum()
        self.ne_ave = (self.ne*MESH.isPlasma).sum()/sum_isPlasma
        self.ni_ave = (self.ni*MESH.isPlasma).sum()/sum_isPlasma
        self.Te_ave = (self.Te*MESH.isPlasma).sum()/sum_isPlasma
        self.Ti_ave = (self.Ti*MESH.isPlasma).sum()/sum_isPlasma
        self.pwr_in_ave = (self.pwr_in*MESH.isPlasma).sum()/sum_isPlasma
        self.pwr_in_tot = self.pwr_in_ave * MESH.area
    
    def _calc_coll(self):
        """Calc collision frequency."""
        temp_coll =self._press/10.0*1e7
        self.coll_em = np.ones_like(self.ne)*temp_coll
        self.coll_im = np.ones_like(self.ne)*temp_coll
    
    def get_eps(self, MESH):
        """Get epsilon from materials."""
        self.eps = np.ones_like(self.ne)
        for idx, mat in np.ndenumerate(MESH.mat):
            if mat == 2:
                self.eps[idx] = 3.8
                
    def savez(self, fname):
        """Save attributes to a bin file."""
        np.savez(fname, 
                 self._press, self._Me, self._Mi, self._wrf,
                 self.ne, self.ni, self.nn,
                 self.Te, self.Ti, self.Tn,
                 self.pot)
    