import numpy as np

####################################################################
############### Lambert/Cosine Anglular distribution ###############
###################### f(theta) = cos(theta) #######################
####################################################################

####################################################################
####################### Cosine SCIPY ###############################
####################################################################

from scipy.stats import cosine
theta_scipy = cosine.rvs(size=10000)
theta_scipy = theta_scipy/2.
theta_scipy = np.degrees(theta_scipy)
# theta_scipy = np.abs(theta_scipy)

####################################################################
####################### Cosine SCIPY ###############################
####################################################################

from scipy.stats import uniform
rv = uniform.rvs(size=10000)
sin_theta = np.sqrt( rv )
# cos_theta = np.sqrt(1.0 - sin_theta * sin_theta)
cos_theta = np.sqrt( 1.0 - rv )
# theta_sqrt = np.arccos( cos_theta )
theta_sqrt = np.arccos( rv )
# theta_sqrt = np.degrees(theta_sqrt)

####################################################################
####################### Cosine Self-Defined ########################
####################################################################

from scipy.stats import rv_continuous
class cosine_func(rv_continuous):
    "Maxwell2d distribution"
    def _pdf(self, x):
        return np.cos(x)

from packages.Constants import PI
cosine_self = cosine_func(momtype=1, a=-PI/2., b=PI/2.)
theta_self = cosine_self.rvs(size=10000)
theta_self = np.degrees(theta_self)
theta_self = np.abs(theta_self)

####################################################################
####################### Cosine Plot ################################
####################################################################

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)
ax.hist(theta_scipy, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(theta_sqrt, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.hist(theta_self, bins=100, density=True, 
        histtype='stepfilled', facecolor='blue', 
        alpha=0.8)
ax.legend(['scipy', 'sqrt', 'self'])
fig.suptitle('Velocity Magitude Distribution in 3D')
plt.show()

fig, ax = plt.subplots(1, 1)
ax.hist(rv, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(cos_theta, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.hist(sin_theta, bins=100, density=True, 
        histtype='step', edgecolor='blue', linewidth=2,
        alpha=0.8)
ax.legend(['rv', 'cos', 'sin'])
fig.suptitle('Velocity Magitude Distribution in 3D')
plt.show()

fig, ax = plt.subplots(1, 1)
ax.hist(theta_sqrt, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(theta_sqrt, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.hist(theta_sqrt, bins=100, density=True, 
        histtype='step', edgecolor='blue', linewidth=2,
        alpha=0.8)
ax.legend(['rv', 'cos', 'sin'])
fig.suptitle('Velocity Magitude Distribution in 3D')
plt.show()