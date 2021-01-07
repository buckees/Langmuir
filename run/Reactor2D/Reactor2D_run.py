"""
Plasma model run file.

example.
"""

import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

import matplotlib.pyplot as plt
from copy import deepcopy

########## import Langmuir modules ##########
from packages.Model.Common.Parameter import PARAMETER
from packages.Model.Reactor2D.Reactor2D_mesh import MESH2D
from packages.Model.Reactor2D.Reactor2D_plasma import PLASMA2D
from packages.Model.Reactor2D.Reactor2D_transp import AMBI2D
from packages.Model.Reactor2D.Reactor2D_rct import REACT2D
from packages.Model.Reactor2D.Reactor2D_eergy import EERGY2D
from packages.Model.Reactor2D.Reactor2D_field import FIELD2D
from packages.Model.Reactor2D.Reactor2D_main import MAIN

# init operation parameters
oper = PARAMETER()
oper.num_iter = 300
oper.dt = 1e-6
oper.ne = 1e18
oper.Te = 2.0

# init mesh obj
fname = 'ICP2D_Mesh'
mesh = MESH2D(fname)
mesh.readin_mesh(fname)

# init plasma obj
pla = PLASMA2D('Plasma')

# init transport module
txp = AMBI2D('Ambipolar')

MAIN(oper, mesh, pla, txp)