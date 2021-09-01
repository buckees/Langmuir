# -*- coding: utf-8 -*-
"""
2D Reactor Model
"""

import numpy as np

eCharge = 1.602e-19 # in C
eMass = 9.109e-31 # in kg
mass_Ar = 6.64e-26 # in kg
eps = 8.85419e-12 # vacuum permitivitty

width = 0.4 # in m
height = 0.4 # in m
depth = 0.4 # in m

res_x = 0.02 # in m
res_z = 0.02 # in m

dt_ne = 1e-5 # in s
dt_Te = 1e-5 # in s

scale = 1

n_iter = 1e4
n_sp = 2
sp_charge = [-1, 1]

freq_coll = 2e7 # in Hz
freq_rf = 13.56e6 # in Hz

vol_plasma = width * height * depth # in m3

den_Ar = 6.64e20 # in m3, 10 mTorr
den_e = 1e15 # in m3
den_Arp = 1e15 # in m3
den_init = [den_e, den_Arp]

mass = [eMass, mass_Ar]

mu=[]
for i in range(n_sp):
    mu.append( eCharge * sp_charge[i] / mass[i] / freq_coll )
    
kT_e = 1.0 # in eV
kT_Arp = 0.1 # in eV
kT = [kT_e, kT_Arp]