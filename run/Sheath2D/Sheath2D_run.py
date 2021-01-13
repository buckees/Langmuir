"""Sheath Model 2D. Main program."""

import os
import glob

import numpy as np
from math import sin
import matplotlib.pyplot as plt

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Sheath2D.Sheath2D_main import MAIN
from packages.Model.Sheath2D.Sheath2D_field import FIELD_SHEATH
from packages.Model.Common.Particle_Mover import EULER_MOVE, LEAPFROG
from packages.Constants import PI

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 10000
oper.max_step = 1000
oper.d_sh = 0.01
oper.wfr_loc = 0.0
oper.dt = 1e-8
oper.Vdc = 100.0
oper.Vrf = 50.0
oper.freq = 2e6
oper.iplot = False

def Efunc_single(t):
    d_sh = oper.d_sh
    freq = oper.freq
    Vdc = oper.Vdc
    Vrf = oper.Vrf
    return -Vdc/d_sh - Vrf/d_sh*sin(2*PI*freq*t)
    

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

field = FIELD_SHEATH('Sheath')
field.add_Efunc(Efunc_single)

erg, ang = MAIN(oper, ptcl, field, move=LEAPFROG)

fname = f'freq{int(oper.freq/1e6)}_Vdc{int(oper.Vdc)}_Vrf{int(oper.Vrf)}'

for i in glob.glob(fname + '.*'):
    os.remove(i)

np.save(fname, erg)

fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                           constrained_layout=True)
ax.hist(erg, bins=50, density=True)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
fig.savefig(fname + '.png', dpi=600)
plt.close()

