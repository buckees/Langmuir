"""Sheath Model 2D. Main program."""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, sqrt

fpng = 'distribution'
vel = np.load('vel.npy')
erg = np.load('erg.npy')
ang = np.load('ang.npy')

vel2d = np.delete(vel, 2, 1)
np.save('vel2d.npy', vel2d)

fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=600,
                           constrained_layout=True)

ax = axes[0, 0]
ax.hist2d(vel[:, 0], vel[:, 1], bins=200, density=False)
ax.set_title('Velocity Distribution')
ax.set_xlabel('Vx (m/s)')
ax.set_ylabel('Vz (m/s)')

ax = axes[0, 1]
ax.hist2d(ang, erg, bins=200, density=False)
ax.set_title('Velocity Distribution')
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Energy (eV)')

ax = axes[1, 0]
ax.hist(erg, bins=100, density=False)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
ax.set_xlim([0, 200])

ax = axes[1, 1]
ax.hist(ang, bins=100, density=False)
ax.set_title('Ion Angular Distribution')
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Count')
ax.set_xlim([-4, 4])


fig.savefig(fpng + '.png', dpi=600)
plt.close()

