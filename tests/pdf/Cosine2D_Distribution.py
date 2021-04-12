import numpy as np
from math import pi, cos

####################################################################
############### Lambert/Cosine Anglular distribution ###############
###################### f(theta) = cos(theta) #######################
####################################################################

size = 10000

####################################################################
####################### Cosine SCIPY ###############################
####################################################################

from scipy.stats import cosine
theta_scipy = cosine.rvs(size=size)
theta_scipy = theta_scipy
theta_scipy = np.degrees(theta_scipy)

####################################################################
####################### Cosine INVERSE #############################
####################################################################

from scipy.stats import uniform
rv = uniform.rvs(size=size)
rv = rv * 2.0 -1.0
theta_inv = np.arccos( rv )
theta_inv = np.degrees(theta_inv)

####################################################################
####################### Cosine Self-Defined ########################
####################################################################

from scipy.stats import rv_continuous
class cosine_func(rv_continuous):
    "Maxwell2d distribution"
    def _pdf(self, x):
        return ( 1.0 + np.cos(x) )/( 2 * pi )

cosine_self = cosine_func(momtype=1, a=-pi, b=pi)
theta_self = cosine_self.rvs(size=size)
theta_self = np.degrees(theta_self)

####################################################################
####################### Cosine MC ##################################
####################################################################

import random
theta_mc = list()
num_count = 0
while num_count < size:
    rv1 = random.uniform(-pi, pi)
    rv2 = random.uniform(0.0, 2.0)
    thres = 1.0 + cos( rv1 )
    if rv2 <= thres:
        theta_mc.append(rv1)
        num_count += 1
theta_mc = np.degrees(theta_mc)

####################################################################
####################### Cosine Plot ################################
####################################################################

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)
ax.hist(theta_scipy, bins=100, density=True, 
        histtype='step', edgecolor='red', linewidth=2,
        alpha=0.8)
ax.hist(theta_inv, bins=100, density=True, 
        histtype='step', edgecolor='green', linewidth=2,
        alpha=0.8)
ax.hist(theta_self, bins=100, density=True, 
        histtype='step', edgecolor='blue', linewidth=2,
        alpha=0.8)
ax.hist(theta_mc, bins=100, density=True, 
        histtype='step', edgecolor='black', linewidth=2,
        alpha=0.8)
ax.legend(['scipy', 'inv', 'self', 'mc'])
fig.suptitle('Cosine Distribution in 2D')
plt.show()
