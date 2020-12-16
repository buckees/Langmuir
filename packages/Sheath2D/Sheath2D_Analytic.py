"""Sheath Model 2D. Main program."""

import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from packages.Constants import (PI, UNIT_CHARGE, AMU)

# Vsh(t) = Vdc + Vrf*sin(wt)
Vdc = 10  # V, average dc voltage
Vrf = 10  # V, average rf voltage amplitude
Lsh = 10e-3  # m, sheath thickness
freq = 13.56e6  # Hz, RF frequency
w = 2.0*PI*freq  # rads/s, RF angular frequency
M = 40  # AMU, Argon mass

dEi = 2.0*UNIT_CHARGE*Vrf/w/Lsh
dEi *= sqrt(2.0*UNIT_CHARGE*Vdc/(M*AMU))


E = np.linspace(0, 1e3, 1001)

prob = 2.0/w/dEi
temp = 4.0/dEi**2
temp *= np.power(E - UNIT_CHARGE*Vdc, 2)
prob /= np.sqrt(1 - temp)


    