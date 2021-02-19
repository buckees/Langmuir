"""Sheath Model 2D. Main program."""

import numpy as np
from copy import copy, deepcopy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
colMap = copy(cm.get_cmap("jet"))
colMap.set_under(color='white')

def ANALYSIS(oper, mesh, pla):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    pla: PLASMA2D(obj), contains all plasma parameters.
    transp: TRANSP2D(obj), contains all transport information.
    eergy: EERGY2D(obj), contains all collision information.
    """

    ########## init and plot plasma ##########
    if oper.irestart:
        pla.loadz(oper.frestart)
        pla.update_plasma(mesh)
        print('restart file is loaded.')
    else:
        return print('oper.irestart is not True, program aborted !')
        
    mesh.plot_var(var=[pla.ne, pla.ni], 
                  var_name=['E Density', 'Ion Density'],
                  fname='Density.png')
    mesh.plot_var(var=[pla.Te, pla.Ti], 
                  var_name=['E Temperature', 'Ion Temperature'],
                  fname='Temperature.png')
    
    dnex, dnez = mesh.cnt_diff(pla.ne)
    dnix, dniz = mesh.cnt_diff(pla.ni)
    mesh.plot_var(var=[dnex, dnez], 
                  var_name=['Delta E Density', 'Delta Ion Density'],
                  fname='Delta_Density.png')
    
    Ex, Ez = dnex/pla.ne, dnez/pla.ne
    for idx, mat in np.ndenumerate(mesh.mat):
           if mat:
               Ex[idx], Ez[idx] = 0.0, 0.0
    Etot = np.sqrt(Ex**2 + Ez**2)
    mesh.plot_var(var=[Ex, Ez, Etot], 
                  var_name=['Ex', 'Ez', 'Etot'],
                  figsize=(24, 8),
                  fname='E-field.png')
    
    print('ANALYSIS program has completed!')
    return Ex, Ez, Etot

########## import Langmuir modules ##########
from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Reactor2D.Reactor2D_mesh import MESH2D
from packages.Model.Reactor2D.Reactor2D_plasma import PLASMA2D
oper = PARAMETER()
oper.irestart = True
oper.frestart = 'ss_10mT_75W.npz'
# init mesh obj
fname = 'ICP2D_Mesh'
mesh = MESH2D(fname)
mesh.readin_mesh(fname)

# init plasma obj
pla = PLASMA2D('Plasma') 

Ex, Ez, Etot = ANALYSIS(oper, mesh, pla)