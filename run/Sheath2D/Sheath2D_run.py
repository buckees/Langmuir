"""Sheath Model 2D. Main program."""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Sheath2D.Sheath2D_main import MAIN
from packages.Model.Sheath2D.Sheath2D_field import FIELD_SHEATH
from packages.Model.Common.Particle_Mover import EULER_MOVE, LEAPFROG

from Efunc import EFUNC

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 2000
oper.max_step = 1000
oper.d_sh = 0.002
oper.wfr_loc = 0.0
oper.imode_move = 'LEAPFROG'
oper.dt = 1e-9
oper.imode_Efunc = 'Dual'
oper.Vdc = 100.0
oper.Vrf = 50.0
oper.freq = 2e6
oper.Vrf2 = 25.0
oper.freq2 = 14e6
oper.iplot = False


# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

field = FIELD_SHEATH('Sheath')
Efunc = EFUNC()
Efunc.load(oper)
if oper.imode_Efunc == "Single":
    field.add_Efunc(Efunc.sgl_freq)
elif oper.imode_Efunc == "Dual":
    field.add_Efunc(Efunc.dual_freq)
elif oper.imode_Efunc == "Customize":
    field.add_Efunc(Efunc.customize)

if oper.imode_move == 'EULER':
    move = EULER_MOVE
elif oper.imode_move == 'LEAPFROG':
    move = LEAPFROG

erg, ang = MAIN(oper, ptcl, field, move=move)

fname = 'test'
# fname = f'freq{int(oper.freq/1e6)}_Vdc{int(oper.Vdc)}_Vrf{int(oper.Vrf)}'
# fname += '_H2O'
# fname = f'dual_freq{int(oper.freq/1e6)}_Vdc{int(oper.Vdc)}_Vrf{int(oper.Vrf)}'
# fname = 'dual_freq214_Vdc100_Vrf5010'

for i in glob.glob(fname + '.*'):
    os.remove(i)

np.save(fname, erg)

fig, axes = plt.subplots(1, 2, figsize=(8, 3), dpi=600,
                           constrained_layout=True)

ax = axes[0]
ax.hist(erg, bins=100, density=False)
ax.set_title('Ion Energy Distribution')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Count')

ax = axes[1]
ax.hist(ang, bins=100, density=False)
ax.set_title('Ion Angular Distribution')
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Count')

fig.savefig(fname + '.png', dpi=600)
plt.close()

