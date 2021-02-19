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
    
    dne = mesh.cnt_diff(self.ne)
    dni = mesh.cnt_diff(self.ni)
    mesh.plot_var(var=[pla.ne, pla.ni], 
                  var_name=['Delta E Density', 'Delta Ion Density'],
                  fname='Delta_Density.png')
    
    
    print('ANALYSIS program has completed!')