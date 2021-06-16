####################################################################
######################### Vector 3D Sampling #######################
####################################################################

import numpy as np

####################################################################
###################### Cartesian Coordinate ########################
####################################################################

from scipy.stats import uniform

uniform_3d = uniform(loc=-1, scale=2)
vec_x = uniform_3d.rvs(size=10000)
vec_z = uniform_3d.rvs(size=10000)
vec_y = uniform_3d.rvs(size=10000)

vec_mag = np.sqrt( vec_x * vec_x + vec_z * vec_z + vec_y * vec_y )
vec_x = np.divide(vec_x, vec_mag, 
                  out=np.zeros_like(vec_mag), where=vec_mag!=0.0)
vec_z = np.divide(vec_z, vec_mag, 
                  out=np.zeros_like(vec_mag), where=vec_mag!=0.0)
vec_y = np.divide(vec_y, vec_mag, 
                  out=np.zeros_like(vec_mag), where=vec_mag!=0.0)


####################################################################
###################### Plot ########################
####################################################################

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)
ax.hist(vec_x, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(vec_z, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.hist(vec_y, bins=100, density=True, 
        histtype='step', edgecolor='blue', linewidth=2,
        alpha=0.8)
ax.legend(['vec_x', 'vec_z', 'vec_y'])
fig.suptitle('Velocity Direction Distribution in 3D')
plt.show()