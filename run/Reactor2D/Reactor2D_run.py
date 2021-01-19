"""
Plasma model run file.

example.
"""

import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

import numpy as np

########## import Langmuir modules ##########
from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Reactor2D.Reactor2D_mesh import MESH2D
from packages.Model.Reactor2D.Reactor2D_plasma import PLASMA2D
from packages.Model.Reactor2D.Reactor2D_transp import AMBI2D
from packages.Model.Reactor2D.Reactor2D_rct import REACT2D
from packages.Model.Reactor2D.Reactor2D_eergy import EERGY2D
from packages.Model.Reactor2D.Reactor2D_field import FIELD2D
from packages.Model.Reactor2D.Reactor2D_main import MAIN

# init operation parameters
oper = PARAMETER()
oper.press = 20
oper.num_iter = 50000
oper.num_plot = 50 
oper.dt = 2e-6
oper.ne = 1e19
oper.num_iter_Te = 10
oper.Te = 2.0
oper.input_pwr = 100.0
oper.idiag = True
oper.irestart = False
oper.frestart = 'restart.npz'

# init mesh obj
fname = 'ICP2D_Mesh'
mesh = MESH2D(fname)
mesh.readin_mesh(fname)

# init plasma obj
pla = PLASMA2D('Plasma')

# init transport module
txp = AMBI2D('Ambipolar')

# init eon energy module
eergy = EERGY2D('Eon_Energy')

# init field module
field = FIELD2D('Field')

def Ey_func(x, z):
    """
    Produce E-field in y direction as a function of space.
    
    Ey(x, z) = Ey_const * (exp(-dist((x, z) to (x1, z1) )/decay_const) + 
                           exp(-dist((x, z) to (x2, z2) )/decay_const))
    x, z: 2d matrix 
    Ey: 2d matrix as x and z, E-field in y direction.
    """
    Ey_const = 7e-1
    decay_const = 0.2
    x1, z1 = -0.12, 0.25
    x2, z2 =  0.12, 0.25
    dist1 = np.sqrt((x - x1)**2 + (z - z1)**2)
    dist2 = np.sqrt((x - x2)**2 + (z - z2)**2)
    Ey = Ey_const * (np.exp(-dist1/decay_const) + 
                     np.exp(-dist2/decay_const))
    return Ey

field.add_Efunc(Ey_func, mesh)

# mesh.plot_var(var=[field.Ey, field.Ey], 
#                   var_name=['Ey', 'Ey'],
#                   fname='init_Ey.png')

# init reaction module
rct = REACT2D('React')

MAIN(oper, mesh, pla, txp, eergy, rct, field)