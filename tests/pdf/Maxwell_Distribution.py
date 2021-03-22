"""
Unit Test for PDF
"""
import numpy as np 
import matplotlib.pyplot as plt

from packages.Constants import PI, EV2J, AMU, K2J, EON_MASS


from scipy.stats import maxwell
from numpy.random import normal, uniform, randint


r = maxwell.rvs(size=1000)

x = np.linspace(maxwell.ppf(0.01),
                maxwell.ppf(0.99), 100)

fig, ax = plt.subplots(1, 1)
ax.plot(x, maxwell.pdf(x),
       'r-', lw=5, alpha=0.6, label='maxwell pdf')

ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.show()



Te = 2.0 # eon temperature in eV
kTe = Te * EV2J
mass_e = EON_MASS

sigma = kTe/mass_e

#################### MAXWELL 3D ####################################

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
fig.suptitle('Velocity Magitude Distribution in 3D')
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
ax.set_title('Velocity in X direction')

ax = axes[1]
ax.hist(vely_maxwell3d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(vely, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)
ax.set_title('Velocity in Z direction')

ax = axes[2]
ax.hist(velz_maxwell3d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(velz, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)
ax.set_title('Velocity in Y direction')

plt.show()

#################### MAXWELL 2D ####################################
vel_phi = 2*randint(0, 2, 10000) - 1
velx_maxwell2d = vel_maxwell3d * np.sin(vel_theta)*vel_phi
velz_maxwell2d = vel_maxwell3d * np.cos(vel_theta)

fig, axes = plt.subplots(1, 2, figsize=(8, 6), dpi=600,
                                constrained_layout=True)

plt.setp(axes, xlim=[-2e12, 2e12], ylim=[0.0, 3e-12])

ax = axes[0]
ax.hist(velx_maxwell2d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(velx, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)
ax.set_title('Velocity in X direction')

ax = axes[1]
ax.hist(velz_maxwell2d, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue',
        alpha=0.5)
ax.hist(vely, bins=100, density=True, 
        histtype='step', edgecolor='red', 
        alpha=0.5)
ax.set_title('Velocity in Z direction')

plt.show()
