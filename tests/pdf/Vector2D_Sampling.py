####################################################################
######################### Vector 2D Sampling #######################
####################################################################

import numpy as np
from math import pi
from scipy.stats import uniform
import matplotlib.pyplot as plt

####################################################################
###################### Cartesian Coordinate ########################
####################################################################

# This cartesian is actually exactly the same as polar.
uniform_cartesian = uniform( loc=-pi, scale=2*pi )
rv = uniform_cartesian.rvs( size=10000 )
x_cart = np.cos( rv )
z_cart = np.sqrt( 1.0 - x_cart * x_cart )
phi_cart = np.arccos( x_cart )

####################################################################
###################### Plot Cartesian ##############################
####################################################################

fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=600,
                         constrained_layout=True)

ax = axes[0]
ax.hist(x_cart, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(z_cart, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.legend(['x_cart', 'z_cart'])
ax.set_title('Velocity Direction Distribution in 2D')

ax = axes[1]
ax.hist(phi_cart, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.set_title('Angle Distribution in 2D')
plt.show()

# fig, ax = plt.subplots(1, 1)
# ax.scatter(vec_x, vec_z)
# plt.show()

####################################################################
########################## Polar Coordinate ########################
####################################################################

uniform_polar = uniform(loc=-pi, scale=2*pi)
phi_polar = uniform_polar.rvs(size=10000)
x_polar = np.cos( phi_polar )
z_polar = np.sin( phi_polar )

####################################################################
########################## Plot Polar ##############################
####################################################################

fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=600,
                         constrained_layout=True)
ax = axes[0]
ax.hist(x_polar, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(z_polar, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.legend(['x_polar', 'z_polar'])
ax.set_title('Velocity Direction Distribution in 2D')

ax = axes[1]
ax.hist(np.degrees(phi_polar), bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.set_title('Angle Distribution in 2D')
plt.show()

# fig, ax = plt.subplots(1, 1)
# ax.scatter(x_polar, z_polar)
# plt.show()

####################################################################
###################### Cartesian Coordinate WRONG ##################
####################################################################

uniform_cartesian = uniform(loc=-1, scale=2)
vec_x = uniform_cartesian.rvs(size=10000)
vec_z = uniform_cartesian.rvs(size=10000)

vec_mag = np.sqrt( vec_x * vec_x + vec_z * vec_z )
vec_x = np.divide(vec_x, vec_mag, 
                  out=np.zeros_like(vec_mag), where=vec_mag!=0.0)
vec_z = np.divide(vec_z, vec_mag, 
                  out=np.zeros_like(vec_mag), where=vec_mag!=0.0)
vec_phi = np.arccos( vec_z )

####################################################################
###################### Plot Cartesian WRONG ########################
####################################################################

fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=600,
                         constrained_layout=True)

ax = axes[0]
ax.hist(vec_x, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(vec_z, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.legend(['vec_x', 'vec_z'])
ax.set_title('Velocity Direction Distribution in 2D')

ax = axes[1]
ax.hist(vec_phi, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.set_title('Angle Distribution in 2D')
plt.show()

# fig, ax = plt.subplots(1, 1)
# ax.scatter(vec_x, vec_z)
# plt.show()
