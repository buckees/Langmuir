"""Sheath Model 2D. Main program."""

import os
import glob

import numpy as np
from math import sin
import matplotlib.pyplot as plt

fname_2 = 'freq2_Vdc100_Vrf50.npy'
fname_14 = 'freq14_Vdc100_Vrf50.npy'

erg2 = np.load(fname_2)
erg14 = np.load(fname_14)


kwargs = dict(histtype='stepfilled', alpha=0.5, density=True, bins=100, ec="k")
fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                           constrained_layout=True)
ax.hist(erg2, **kwargs)
ax.hist(erg14, **kwargs)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Density')
fig.savefig('freq' + '.png', dpi=600)
plt.close()

fname_H3 = 'freq14_Vdc100_Vrf50_H3.npy'
fname_H2O = 'freq14_Vdc100_Vrf50_H2O.npy'
fname_Ar = 'freq14_Vdc100_Vrf50.npy'

ergH3 = np.load(fname_H3)
ergH2O = np.load(fname_H2O)
ergAr = np.load(fname_Ar)

fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                           constrained_layout=True)
ax.hist(ergH3, **kwargs)
ax.hist(ergH2O, **kwargs)
ax.hist(ergAr, **kwargs)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Density')
fig.savefig('mass' + '.png', dpi=600)
plt.close()


fname_dual1 = 'dual_freq214_Vdc100_Vrf5010.npy'
fname_dual2 = 'dual_freq214_Vdc100_Vrf5050.npy'

ergDual1 = np.load(fname_dual1)
ergDual2 = np.load(fname_dual2)

fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                           constrained_layout=True)
ax.hist(erg2, **kwargs)
ax.hist(ergDual1, **kwargs)
ax.hist(ergDual2, **kwargs)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Density')
fig.savefig('dual' + '.png', dpi=600)
plt.close()
