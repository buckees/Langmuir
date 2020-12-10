"""
Plasma model run file.

example.
"""

import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

# import numpy as np


"""Import plasma modules."""
from packages.Reactor2D.Reactor2D_Mesh import MESH2D
from packages.Reactor2D.Reactor2D_Plasma import PLASMA2D
from packages.Reactor2D.Reactor2D_Transp import AMBI2D
from packages.Reactor2D.Reactor2D_React import REACT2D
from packages.Reactor2D.Reactor2D_Eergy import EERGY2D

fname = 'ICP2D_Mesh'
# init MESHGRID obj
MESH = MESH2D(fname)
# readin mesh
MESH.readin_mesh(fname)

PLA = PLASMA2D('PlA2D')
PLA.import_mesh(MESH)
PLA.init_plasma(ne=1e16, Te=1.5)

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

'Here!'

# init TRANSP2D module
TXP = AMBI2D('TXP2D')
TXP.from_PLASMA(PLA)
TXP.calc_ambi(PLA)

# init React module
SRC = REACT2D('SRC2D')
SRC.from_PLASMA(PLA)
SRC.calc_src(PLA)

# init Eergy module
EERN = EERGY2D('EERN2D')
EERN.from_PLASMA(PLA)
EERN.calc_Te(PLA)

dt = 1e-5
niter = 100
for itn in range(niter):
    TXP.calc_ambi(PLA)
    TXP.solve_fluid(dt)
    if not (itn+1) % (niter/10):
        TXP.to_PLASMA(PLA)
        MESH.plot_var(var=[PLA.ne, PLA.ni], 
              var_name=['E Density', 'Ion Density'],
              fname=f'plasma_itn{itn+1}')
        TXP.from_PLASMA(PLA)
        
        

# ne_ave, ni_ave, Te_ave = [], [], []
# time = []
# niter = 100
# dt = 1e-5
# niter_Te = 30
# for itn in range(niter):
#     pwr2d.calc_pwr_in(pla2d, pwr=100.0, imode='ne')
#     txp2d.calc_ambi(pla2d)
#     een2d.get_pwr(pwr2d)
#     for itn_Te in range(niter_Te):    
#         een2d.calc_Te(dt/niter_Te, pla2d, txp2d)
#     pla2d.get_Te(een2d)
#     src2d.calc_src(pla2d)
#     pla2d.den_evolve(dt, txp2d, src2d)
#     ne_ave.append(pla2d.ne.mean())
#     ni_ave.append(pla2d.ni.mean())
#     Te_ave.append(pla2d.Te.mean())
#     time.append(dt*(itn+1))
#     if not (itn+1) % (niter/10):
#         pla2d.plot_plasma(fname=f'plasma_itn{itn+1}', 
#                           figsize=figsize, ihoriz=ihoriz)
#         pla2d.plot_Te(fname=f'Te_itn{itn+1}', 
#                           figsize=figsize, ihoriz=ihoriz)
#         een2d.plot_dQe(pla=pla2d, fname=f'dQe_itn{itn+1}', 
#                           figsize=figsize, ihoriz=ihoriz)
#         src2d.plot_src(pla=pla2d, fname=f'src_itn{itn+1}', 
#                           figsize=figsize, ihoriz=ihoriz)
#         txp2d.plot_flux(pla=pla2d, fname=f'flux_itn{itn+1}',
#                         figsize=figsize, ihoriz=ihoriz)
#         pla2d.calc_conde(2*PI*13.56e6)
#         pla2d.plot_conde(fname=f'conde_itn{itn+1}', 
#                           figsize=figsize, ihoriz=ihoriz)

#         # plot ave. values
#         fig, axes = plt.subplots(1, 2, figsize=(8,4), dpi=300,
#                                              constrained_layout=True)
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
