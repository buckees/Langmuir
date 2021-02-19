"""Feature Model 2D. Main program."""

import numpy as np
import random
import matplotlib.pyplot as plt

from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Common.Particle import PARTICLE
from packages.Model.Feature2D.Feature2D_mesh import MESH2D
from packages.Model.Feature2D.Feature2D_rflct import REFLECT
from packages.Model.Feature2D.Feature2D_rct import REACT
from packages.Model.Feature2D.Feature2D_main import MAIN

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 10000
oper.max_step = 1000
oper.step_fac = 0.5
oper.max_rflct = 5
oper.surf_norm_range = 3
oper.surf_norm_mode = 'Sum Vector'
oper.num_plot = 5
oper.prob_rflct = 0.5
oper.idiag = True


# init mesh
mesh = MESH2D()
# readin mesh
fname = 'SiEtch_Base_Mesh'
mesh.readin_mesh(fname)

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

def xFunc():
    return mesh.init_ptcl_posn()

vel_list = np.load('vel.npy')
num = len(vel_list)
def vFunc():
    idx = random.randint(0, num-1)
    return vel_list[idx]

ptcl.add_initPosnFunc(xFunc)
ptcl.add_initVelFunc(vFunc)

# init reflection
rflct = REFLECT()

# init reaction
rct = REACT()

MAIN(oper, ptcl, mesh, rct, rflct)

if oper.idiag:
    vel = np.load('Feat2D_initVel.npy')
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=600,
                               constrained_layout=True)
    
    ax = axes[0, 0]
    ax.hist2d(vel[:, 0], vel[:, 1], bins=200, density=False)
    ax.set_title('Velocity Distribution')
    ax.set_xlabel('Vx (m/s)')
    ax.set_ylabel('Vz (m/s)')
    
    ax = axes[0, 1]
    ax.hist2d(vel[:, 0], vel[:, 1], bins=200, density=False)
    ax.set_title('Velocity Distribution')
    ax.set_xlabel('Angle (degree)')
    ax.set_ylabel('Energy (eV)')
    
    ax = axes[1, 0]
    ax.hist(vel[:, 0], bins=100, density=False)
    ax.set_title('Ion Energy Distribution')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Count')
    # ax.set_xlim([0, 200])
    
    ax = axes[1, 1]
    ax.hist(vel[:, 1], bins=100, density=False)
    ax.set_title('Ion Angular Distribution')
    ax.set_xlabel('Angle (degree)')
    ax.set_ylabel('Count')
    # ax.set_xlim([-4, 4])
    
    fpng = 'init_distrb'
    fig.savefig(fpng + '.png', dpi=600)
    plt.close()
    
