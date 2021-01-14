"""Sheath Model 2D. Main program."""

import numpy as np
from copy import copy, deepcopy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
colMap = copy(cm.get_cmap("jet"))
colMap.set_under(color='white')

def MAIN(oper, mesh, pla, txp, eergy=None, rct=None, field=None):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    pla: PLASMA2D(obj), contains all plasma parameters.
    transp: TRANSP2D(obj), contains all transport information.
    eergy: EERGY2D(obj), contains all collision information.
    """

    ########## init and plot plasma ##########
    pla.init_plasma(MESH=mesh, ne=oper.ne, Te=oper.Te)
    if oper.idiag:        
        mesh.plot_var(var=[pla.ne, pla.ni], 
                      var_name=['E Density', 'Ion Density'],
                      fname='Init_Density.png')
        mesh.plot_var(var=[pla.Te, pla.Ti], 
                      var_name=['E Temperature', 'Ion Temperature'],
                      fname='Init_Temperature.png')
    ##########################################

    dt = oper.dt
    
    ########## pre-run tranport ##########
    txp.from_PLASMA(pla)
    
    for itn in range(500):
        txp.calc_ambi(mesh)
        txp.solve_fluid(dt)

    txp.to_PLASMA(pla)
    pla.update_plasma(mesh)
    if oper.idiag:    
        mesh.plot_var(var=[pla.ne, pla.ni], 
                      var_name=['E Density', 'Ion Density'],
                      fname='Prerun_Density.png')
    
    ####################################
            
    ########## init and plot field ##########
    field.to_PLASMA(pla)
    
    mesh.plot_var(var=[pla.Ey, pla.Ex], 
          var_name=['Ey', 'Ex'],
          fname='init_EField')
    #########################################
    
    ########## pre-run eon energy ##########
    eergy.from_PLASMA(pla)
    for itn in range(100):
        eergy.solve_Te(mesh, dt)
        
    eergy.to_PLASMA(pla)
    pla.update_plasma(mesh)
    if oper.idiag:    
        mesh.plot_var(var=[pla.Te, pla.Ti], 
                      var_name=['E Temperature', 'Ion Temperature'],
                      fname='Prerun_Temperature.png')
    ########################################
    

    ##########################################
    ########## main loop for plasma ##########
    ##########################################
    ne_ave, ni_ave, Te_ave = [], [], []
    pwr_in_tot = []
    time = []
   
    for itn in range(oper.num_iter):
        
        ########## print progress ##########
        if (itn + 1) % int(oper.num_iter/oper.num_plot) == 0:
            print('%d iterations have completed!' % (itn+1))
        ####################################
        # call field module
        field.from_PLASMA(pla)
        field.adjust_E(oper.input_pwr)
        field.to_PLASMA(pla)
        # call electron energy module
        eergy.from_PLASMA(pla)
        for itn_Te in range(oper.num_iter_Te):    
            eergy.solve_Te(mesh, dt/oper.num_iter_Te)
        eergy.to_PLASMA(pla)
        pla.update_plasma(mesh)
        # call reaction module
        rct.from_PLASMA(pla)
        rct.calc_src(mesh)
        rct.to_PLASMA(pla)
        pla.update_plasma(mesh)
        # call transport module
        txp.from_PLASMA(pla)
        txp.calc_ambi(mesh)
        txp.solve_fluid(dt)
        txp.to_PLASMA(pla)
        pla.update_plasma(mesh)
        
        
        ########## record and plot ##########
        ne_ave.append(deepcopy(pla.ne_ave))
        ni_ave.append(deepcopy(pla.ni_ave))
        Te_ave.append(deepcopy(pla.Te_ave))
        pwr_in_tot.append(deepcopy(pla.pwr_in_tot))
        time.append(dt*(itn+1))
        if oper.idiag:
            if not (itn+1) % (oper.num_iter/oper.num_plot):
                # plot 2D        
                mesh.plot_var(var=[pla.ne, pla.ni], 
                      var_name=['E Density', 'Ion Density'],
                      fname=f'Density_itn{itn+1}')
                mesh.plot_var(var=[pla.Te, pla.Ti], 
                      var_name=['E Temperature', 'Ion Temperature'],
                      fname=f'Te_itn{itn+1}')
                mesh.plot_var(var=[pla.Se, pla.dfluxe], 
                      var_name=['E Prod', 'E Loss'],
                      fname=f'rct_itn{itn+1}')
                mesh.plot_var(var=[pla.pwr_in, eergy.dQe], 
                      var_name=['Power due to Ey', 'dQe'],
                      fname=f'Power_itn{itn+1}')
        
                # plot ave. values
                fig, axes = plt.subplots(1, 3, figsize=(12,4), dpi=300,
                                                      constrained_layout=True)
                ax = axes[0]
                ax.plot(time, ne_ave, 'b-')
                ax.legend(['ne'])
                ax.set_title('Eon Density (m^-3)')
                plt.xlabel('Time (s)')
                plt.ylabel('Ave. Density (m^-3)')
                
                ax = axes[1]
                ax.plot(time, Te_ave, 'r-')
                ax.legend(['Te'])
                ax.set_title('Eon Temperature (eV)')
                plt.xlabel('Time (s)')
                plt.ylabel('Ave. Eon Temperature (eV)')
                
                ax = axes[2]
                ax.plot(time, pwr_in_tot, 'k-')
                ax.legend(['Total Power'])
                ax.set_title('Total Power (W)')
                plt.xlabel('Time (s)')
                plt.ylabel('Total Input Power (W)')
                
                fig.savefig('Ave_vs_Time.png', dpi=300)
                plt.close()
        