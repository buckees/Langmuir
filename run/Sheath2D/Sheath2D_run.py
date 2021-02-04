"""Sheath Model 2D. Main program."""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, sqrt

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Sheath2D.Sheath2D_main import MAIN
from packages.Model.Sheath2D.Sheath2D_field import FIELD_SHEATH
from packages.Model.Common.Particle_Mover import EULER_MOVE, LEAPFROG
from scipy.stats import maxwell
from packages.Constants import EV2J, AMU

from Efunc import EFUNC

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 20000
oper.max_step = 1000
oper.Ti = 1.0  # eV
oper.d_sh = 0.002  # m
oper.wfr_loc = 0.0
oper.imode_move = 'LEAPFROG'
oper.dt = 1e-9
oper.imode_Efunc = 'Dual'
oper.Vdc = 0.0
oper.Vrf = 100.0
oper.freq = 2e6
oper.Vrf2 = 20.0
oper.freq2 = 14e6

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

def xFunc():
    return np.array([0.0, oper.d_sh, 0.0])

def vFunc():
    a = sqrt(oper.Ti*EV2J/(ptcl.mass*AMU))  # a = sqrt(kT/m)
    speed = maxwell.rvs(loc=0.0, scale=a, size=1)
    vel = np.zeros(3)
    mu, sigma = 0.0, 0.1  # default mean and standard deviation
    theta = np.random.normal(mu, sigma)
    vel[0], vel[1] = sin(theta), -cos(theta)
    vel = speed * vel
    return vel

ptcl.add_initPosnFunc(xFunc)
ptcl.add_initVelFunc(vFunc)

# init field
field = FIELD_SHEATH('Sheath')
Efunc = EFUNC()
Efunc.load(oper)
if oper.imode_Efunc == "Single":
    field.add_Efunc(Efunc.sgl_freq)
elif oper.imode_Efunc == "Dual":
    field.add_Efunc(Efunc.dual_freq)
    # freq = min(Efunc.freq, Efunc.freq2)
    # period = 1/freq
    # time = np.linspace(0.0, period, 101)
    # voltage = np.array([Efunc.dual_freq(t) for t in time])
    # plt.plot(time, voltage)
elif oper.imode_Efunc == "Tilt":
    field.add_Efunc(Efunc.tilt)
elif oper.imode_Efunc == "Customize":
    field.add_Efunc(Efunc.customize)

if oper.imode_move == 'EULER':
    move = EULER_MOVE
elif oper.imode_move == 'LEAPFROG':
    move = LEAPFROG

erg, ang, init_erg, init_ang = MAIN(oper, ptcl, field, move=move)

fname = 'test'
# fname = f'freq{int(oper.freq/1e6)}_Vdc{int(oper.Vdc)}_Vrf{int(oper.Vrf)}'
# fname += '_H2O'
# fname = f'dual_freq{int(oper.freq/1e6)}_Vdc{int(oper.Vdc)}_Vrf{int(oper.Vrf)}'
# fname = 'dual_freq214_Vdc100_Vrf5010'

for i in glob.glob(fname + '.*'):
    os.remove(i)

np.save(fname, erg)

fig, axes = plt.subplots(2, 2, figsize=(8, 6), dpi=600,
                           constrained_layout=True)

ax = axes[0, 0]
ax.hist(init_erg, bins=100, density=False)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
ax.set_xlim([0, 12])

ax = axes[0, 1]
ax.hist(init_ang, bins=100, density=False)
ax.set_title('Ion Angular Distribution')
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Count')
ax.set_xlim([-20, 20])

ax = axes[1, 0]
ax.hist(erg, bins=100, density=False)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')
ax.set_xlim([0, 200])

ax = axes[1, 1]
ax.hist(ang, bins=100, density=False)
ax.set_title('Ion Angular Distribution')
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Count')
ax.set_xlim([-4, 4])


fig.savefig(fname + '.png', dpi=600)
plt.close()

