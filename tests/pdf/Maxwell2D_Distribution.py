"""
Unit Test for PDF
"""
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from packages.Constants import PI, EV2J, AMU, K2J, EON_MASS
from numpy.random import normal, uniform, randint

####################################################################
#################### NOT completed!! ###############################
####################################################################


####################################################################
############### Maxwell2D speed distribution #######################
######## f(v) = m/kT * v * exp( -1/2 * m/kT * v**2 ) ###############
####################################################################

Te = 2.0 # eon temperature in eV
kTe = Te * EV2J
mass_e = EON_MASS
sigma = sqrt(kTe/mass_e)

####################################################################
#################### Maxwell2D Self-Defined ########################
####################################################################

from scipy.stats import rv_continuous
class maxwell2d(rv_continuous):
    "Maxwell2d distribution"
    def _pdf(self, x):
        return x * np.exp(- x**2 / 2.)

maxwell_2d = maxwell2d(a=0.0)
vel_self = maxwell_2d.rvs(scale=sigma, size=10000)
phi = uniform(-PI, PI, 10000)
# phi = np.arccos(1.0 - 2.0 * uniform(0.0, 1.0, 10000))

velx_self = vel_self * np.sin(phi)
velz_self = vel_self * np.cos(phi)

####################################################################
#################### Maxwell2D Normal ##############################
####################################################################

velx_norm = normal(scale=sigma, size=10000)
velz_norm = normal(scale=sigma, size=10000)
vel_norm = np.sqrt(velx_norm**2 + velz_norm**2)

####################################################################
#################### Maxwell2D Plot ################################
####################################################################

fig, ax = plt.subplots(1, 1)
ax.hist(vel_self, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(vel_norm, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['self','normal'])
fig.suptitle('Velocity Magitude Distribution in 3D')
plt.show()

fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=600,
                                constrained_layout=True)
# plt.setp(axes, xlim=[-2e12, 2e12], ylim=[0.0, 1.5e-12])

ax = axes[0]
ax.hist(velx_self, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(velx_norm, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['self','normal'])
ax.set_title('Velocity in X direction')

ax = axes[1]
ax.hist(velz_self, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(velz_norm, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['self','normal'])
ax.set_title('Velocity in Z direction')

plt.show()
