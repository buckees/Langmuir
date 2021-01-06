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

"""Import Langmuir modules."""
from packages.Reactor2D.Reactor2D_Mesh import MESH2D
from packages.Reactor2D.Reactor2D_Plasma import PLASMA2D
from packages.Reactor2D.Reactor2D_Transp import AMBI2D
from packages.Reactor2D.Reactor2D_React import REACT2D
from packages.Reactor2D.Reactor2D_Eergy import EERGY2D
from packages.Reactor2D.Reactor2D_Field import FIELD2D

fname = 'ICP2D_Mesh'
# init MESHGRID obj
MESH = MESH2D(fname)
# readin mesh
MESH.readin_mesh(fname)

PLA = PLASMA2D('PlA2D')
PLA.import_mesh(MESH)
PLA.init_plasma(ne=1e18, Te=1.5)

temp_ratio = MESH.width/MESH.height
if temp_ratio > 2.0:
    figsize = (4, 8)
    ihoriz = 0
else:
    figsize = (10, 4)
    ihoriz = 1
    
MESH.plot_var(var=[PLA.ne, PLA.ni], 
              var_name=['E Density', 'Ion Density'],
              fname='Init_Density.png')
MESH.plot_var(var=[PLA.Te, PLA.Ti], 
              var_name=['E Temperature', 'Ion Temperature'],
              fname='Init_Temperature.png')

# init TRANSP2D module
TXP = AMBI2D('TXP2D')
TXP.from_PLASMA(PLA)
TXP.calc_ambi(PLA)


dt = 1e-5
niter = 300
TXP.from_PLASMA(PLA)
for itn in range(niter):
    TXP.calc_ambi(PLA)
    TXP.solve_fluid(dt)
    if not (itn+1) % (niter/5):
        TXP.to_PLASMA(PLA)
        PLA.update_plasma()
        # MESH.plot_var(var=[PLA.ne, PLA.ni], 
        #       var_name=['E Density', 'Ion Density'],
        #       fname=f'Init_itn{itn+1}')
        # MESH.plot_var(var=[PLA.Ex, PLA.Ez], 
        #       var_name=['E-Field in x', 'E-Field in z'],
        #       fname=f'EF_itn{itn+1}')
        TXP.from_PLASMA(PLA)
TXP.from_PLASMA(PLA)

FIELD = FIELD2D('FIELD2D')
FIELD.from_PLASMA(PLA)
FIELD.create_Ey()
FIELD.to_PLASMA(PLA)

MESH.plot_var(var=[PLA.Ey, PLA.Ex], 
      var_name=['Ey', 'Ex'],
      fname='E-Field')

# init Eergy module
EERN = EERGY2D('EERN2D')

dt = 2e-7
niter = 10000
EERN.from_PLASMA(PLA)
for itn in range(niter):
    EERN.solve_Te(PLA, dt)
    if not (itn+1) % (niter/5):
        EERN.to_PLASMA(PLA)
        PLA.update_plasma()
        MESH.plot_var(var=[PLA.Te, PLA.Ti], 
              var_name=['E Temperature', 'Ion Temperature'],
              fname=f'Te_itn{itn+1}')
        MESH.plot_var(var=[PLA.pwr_in, EERN.dQe], 
              var_name=['Power due to Ey', 'dQe'],
              fname=f'Power_itn{itn+1}')
        EERN.from_PLASMA(PLA)
EERN.to_PLASMA(PLA)

# # init React module
# SRC = REACT2D('SRC2D')


# ne_ave, ni_ave, Te_ave = [], [], []  
# time = []
# niter = 1000
# dt = 1e-7
# niter_Te = 40


# for itn in range(niter):
#     # call REACT2D
#     SRC.from_PLASMA(PLA)
#     SRC.calc_src(PLA)
#     SRC.to_PLASMA(PLA)
#     PLA.update_plasma()
#     # call TRANSP2D
#     TXP.from_PLASMA(PLA)
#     TXP.calc_ambi(PLA)
#     TXP.solve_fluid(dt)
#     TXP.to_PLASMA(PLA)
#     PLA.update_plasma()
#     # call EERGY2D
#     EERN.from_PLASMA(PLA)
#     for itn_Te in range(niter_Te):    
#         EERN.solve_Te(PLA, dt/niter_Te)
#     EERN.to_PLASMA(PLA)
#     PLA.update_plasma()
#     # record ave
#     ne_ave.append(deepcopy(PLA.ne_ave))
#     ni_ave.append(deepcopy(PLA.ni_ave))
#     Te_ave.append(deepcopy(PLA.Te_ave))
#     time.append(dt*(itn+1))
#     if not (itn+1) % (niter/5):
#         # plot 2D        
#         MESH.plot_var(var=[PLA.ne, PLA.ni], 
#               var_name=['E Density', 'Ion Density'],
#               fname=f'Density_itn{itn+1}')
#         MESH.plot_var(var=[PLA.Te, PLA.Ti], 
#               var_name=['E Temperature', 'Ion Temperature'],
#               fname=f'Te_itn{itn+1}')
#         MESH.plot_var(var=[PLA.dfluxe, PLA.Se], 
#               var_name=['E Loss', 'E Prod'],
#               fname=f'SRC_itn{itn+1}')
#         MESH.plot_var(var=[PLA.pwr_in, EERN.dQe], 
#               var_name=['Power due to Ey', 'dQe'],
#               fname=f'Power_itn{itn+1}')

#         # plot ave. values
#         fig, axes = plt.subplots(1, 2, figsize=(8,4), dpi=300,
#                                               constrained_layout=True)
#         ax = axes[0]
#         ax.plot(time, ne_ave, 'b-')
#         ax.legend(['ne'])
#         ax.set_title('Eon Density (m^-3)')
#         plt.xlabel('Time (s)')
#         plt.ylabel('Ave. Density (m^-3)')
        
#         ax = axes[1]
#         ax.plot(time, Te_ave, 'r-')
#         ax.legend(['Te'])
#         ax.set_title('Eon Temperature (eV)')
#         plt.xlabel('Time (s)')
#         plt.ylabel('Ave. Eon Temperature (eV)')
        
#         fig.savefig('Ave_vs_Time.png', dpi=300)
#         plt.close()
