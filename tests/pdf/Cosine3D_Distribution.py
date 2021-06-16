"""
Unit Test for PDF
"""
import numpy as np 
from math import sqrt
import matplotlib.pyplot as plt
from packages.Constants import PI, EV2J, AMU, K2J, EON_MASS
from scipy.stats import maxwell
from numpy.random import normal, uniform, randint


####################################################################
############### Maxwell3D speed distribution #######################
# f(theta) = cos(theta)
####################################################################


####################################################################
#################### Cosine3D SCIPY ###############################
####################################################################

vel_scipy = maxwell.rvs(scale=sigma, size=10000)
# vel_theta = uniform(0.0, PI, 10000)
theta = np.arccos(1.0 - 2.0 * uniform(0.0, 1.0, 10000))
phi = uniform(0.0, 2*PI, 10000)

velx_scipy = vel_scipy * np.sin(theta)*np.cos(phi)
vely_scipy = vel_scipy * np.sin(theta)*np.sin(phi)
velz_scipy = vel_scipy * np.cos(theta)

####################################################################
#################### Maxwell3D Normal ##############################
####################################################################

velx_norm = normal(scale=sigma, size=10000)
vely_norm = normal(scale=sigma, size=10000)
velz_norm = normal(scale=sigma, size=10000)
vel_norm = np.sqrt(velx_norm**2 + vely_norm**2 + velz_norm**2)

####################################################################
#################### Maxwell3D Plot ################################
####################################################################

fig, ax = plt.subplots(1, 1)
ax.hist(vel_scipy, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(vel_norm, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['scipy','normal'])
fig.suptitle('Velocity Magitude Distribution in 3D')
plt.show()


fig, axes = plt.subplots(1, 3, figsize=(12 , 4), dpi=600,
                                constrained_layout=True)
# plt.setp(axes, xlim=[-2e12, 2e12], ylim=[0.0, 1.5e-12])

ax = axes[0]
ax.hist(velx_scipy, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(velx_norm, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['scipy','normal'])
ax.set_title('Velocity in X direction')

ax = axes[1]
ax.hist(velz_scipy, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(velz_norm, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['scipy','normal'])
ax.set_title('Velocity in Z direction')

ax = axes[2]
ax.hist(vely_scipy, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(vely_norm, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['scipy','normal'])
ax.set_title('Velocity in Y direction')

plt.show()