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
from packages.Model.Sheath2D.Sheath2D_coll import COLLISION
from packages.Model.Sheath2D.Sheath2D_stats import STATS
from packages.Model.Common.Particle_Mover import EULER_MOVE, LEAPFROG
from scipy.stats import maxwell
from packages.Constants import PI, EV2J, AMU, K2J

from Efunc import EFUNC

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 5000
oper.max_step = 1000
oper.Ti = 1.0  # eV
oper.Tn = 300.0  # K
oper.imode_initVel = 'Cosine'
oper.d_sh = 0.002  # m
oper.wfr_loc = 0.0
oper.imode_move = 'LEAPFROG'
oper.dt = 1e-9
oper.imode_Efunc = 'Single'
oper.Vdc = 100.0
oper.Vrf = 50.0
oper.freq = 14e6
oper.Vrf2 = 20.0
oper.freq2 = 14e6
oper.coll_freq = 4e6
oper.num_report = 5
oper.idiag = True

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

def xFunc():
    return np.array([0.0, oper.d_sh, 0.0])

def vFunc():
    a = sqrt(oper.Ti*EV2J/(ptcl.mass*AMU))  # a = sqrt(kT/m)
    speed = maxwell.rvs(loc=0.0, scale=a, size=1)
    vel = np.zeros(3)
    mu, sigma = 0.0, 10.0  # default mean and standard deviation
    theta = np.random.normal(mu, sigma)
    vel[0], vel[1] = sin(theta), -cos(theta)
    vel = speed * vel
    return vel

ptcl.add_xFunc(xFunc)
ptcl.add_vFunc(vFunc)

# init field
field = FIELD_SHEATH('Sheath')
Efunc = EFUNC()
Efunc.load(oper)
if oper.imode_Efunc == "Single":
    field.add_Efunc(Efunc.sgl_freq)
elif oper.imode_Efunc == "Dual":
    field.add_Efunc(Efunc.dual_freq)
elif oper.imode_Efunc == "Tilt":
    field.add_Efunc(Efunc.tilt)
elif oper.imode_Efunc == "Customize":
    field.add_Efunc(Efunc.customize)

if oper.imode_move == 'EULER':
    move = EULER_MOVE
elif oper.imode_move == 'LEAPFROG':
    move = LEAPFROG
    
def func_CollFreq(ptcl_vel):
    return oper.coll_freq

def func_ReinitVel(ptcl_vel):
    vel = np.zeros(3)
    a = sqrt(oper.Tn*K2J/(ptcl.mass*AMU))  # a = sqrt(kT/m)
    speed = maxwell.rvs(loc=0.0, scale=a, size=1)
    theta = np.random.uniform(-PI, PI)
    vel[0], vel[1] = sin(theta), -cos(theta)
    vel = speed * vel
    return vel

def func_0(ptcl_vel):
    return np.zeros(3)

coll = COLLISION('Ion_Collision')
coll.add_func_CollFreq(func_CollFreq)
coll.add_func_ReinitVel(func_ReinitVel)
# coll.add_func_ReinitVel(func_0)


if oper.idiag:
    stats = STATS()
    vel, stats = MAIN(oper, ptcl, field, coll, move, stats)
else:
    vel = MAIN(oper, ptcl, field, coll, move, stats)

fname = 'IAEDF'
# fname = f'freq{int(oper.freq/1e6)}_Vdc{int(oper.Vdc)}_Vrf{int(oper.Vrf)}'
# fname += '_H2O'
# fname = f'dual_freq{int(oper.freq/1e6)}_Vdc{int(oper.Vdc)}_Vrf{int(oper.Vrf)}'
# fname = 'dual_freq214_Vdc100_Vrf5010'

for i in glob.glob(fname + '.*'):
    os.remove(i)

if oper.idiag:
    stats.save2csv()
    stats.plot()
    
    # fig, axes = plt.subplots(2, 2, figsize=(8, 6), dpi=600,
    #                             constrained_layout=True)
    
    # ax = axes[0, 0]
    # stats.df['Init_Erg'].hist(bins=100, density=False, ax=ax)
    # ax.set_title('Ion Energy Distribution')
    # ax.set_xlabel('Energy (eV)')
    # ax.set_ylabel('Count')
    # # ax.set_xlim([0, 12])
    
    # ax = axes[0, 1]
    # stats.df['Init_Ang'].hist(bins=100, density=False, ax=ax)
    # ax.set_title('Ion Angular Distribution')
    # ax.set_xlabel('Angle (degree)')
    # ax.set_ylabel('Count')
    # # ax.set_xlim([-20, 20])
    
    # ax = axes[1, 0]
    # stats.df['End_Erg'].hist(bins=100, density=False, ax=ax)
    # ax.set_title('Ion Energy Distribution')
    # ax.set_xlabel('Energy (eV)')
    # ax.set_ylabel('Count')
    # # ax.set_xlim([0, 200])
    
    # ax = axes[1, 1]
    # stats.df['End_Ang'].hist(bins=100, density=False, ax=ax)
    # ax.set_title('Ion Angular Distribution')
    # ax.set_xlabel('Angle (degree)')
    # ax.set_ylabel('Count')
    # # ax.set_xlim([-10, 10])
    
    # fig.savefig(fname + '.png', dpi=600)
    # plt.close()
