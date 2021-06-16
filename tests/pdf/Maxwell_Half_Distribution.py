"""
Unit Test for PDF
"""
import numpy as np 
import matplotlib.pyplot as plt

from packages.Constants import PI, EV2J, AMU, K2J, EON_MASS


from scipy.stats import maxwell
from numpy.random import normal, uniform


Te = 2.0 # eon temperature in eV
kTe = Te * EV2J
mass_e = EON_MASS

sigma = kTe/mass_e
vel_maxwell3d = maxwell.rvs(scale=sigma, size=10000)
# vel_theta = uniform(0.0, PI, 10000)
vel_theta = np.arccos(1.0 - 2.0 * uniform(0.0, 1.0, 10000))
vel_phi = uniform(0.0, 2*PI, 10000)

velx_maxwell3d = vel_maxwell3d * np.sin(vel_theta)*np.cos(vel_phi)
vely_maxwell3d = vel_maxwell3d * np.sin(vel_theta)*np.sin(vel_phi)
velz_maxwell3d = vel_maxwell3d * np.cos(vel_theta)

velx = normal(scale=sigma, size=10000)
vely = normal(scale=sigma, size=10000)
velz = normal(scale=sigma, size=10000)

vel = np.sqrt(velx**2 + vely**2 + velz**2)

fig, ax = plt.subplots(1, 1)
ax.hist(vel_maxwell3d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(vel, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)
plt.show()

fig, axes = plt.subplots(1, 3, figsize=(12 , 6), dpi=600,
                                constrained_layout=True)

plt.setp(axes, xlim=[-2e12, 2e12], ylim=[0.0, 3e-12])

ax = axes[0]
ax.hist(velx_maxwell3d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(velx, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)
ax = axes[1]
ax.hist(vely_maxwell3d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(vely, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)
ax = axes[2]
ax.hist(velz_maxwell3d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(velz, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)

plt.show()
