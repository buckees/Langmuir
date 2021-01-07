"""Sheath Model 2D. Main program."""

import numpy as np

def MAIN(oper, mesh, pla, txp, eergy=None, rct=None, field=None):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    pla: PLASMA2D(obj), contains all plasma parameters.
    transp: TRANSP2D(obj), contains all transport information.
    eergy: EERGY2D(obj), contains all collision information.
    """

    ########## init and plot plasma ##########
    pla.import_mesh(mesh)
    pla.init_plasma(ne=oper.ne, Te=oper.Te)
        
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
        txp.calc_ambi(pla)
        txp.solve_fluid(dt)

    txp.to_PLASMA(pla)
    pla.update_plasma()
    
    mesh.plot_var(var=[pla.ne, pla.ni], 
                  var_name=['E Density', 'Ion Density'],
                  fname='Prerun_Density.png')
    
    ####################################
            
    ########## init and plot field ##########
    field.from_PLASMA(pla)
    field.create_Ey()
    field.to_PLASMA(pla)
    
    mesh.plot_var(var=[pla.Ey, pla.Ex], 
          var_name=['Ey', 'Ex'],
          fname='E-Field')
    #########################################
    
    ########## pre-run eon energy ##########
    eergy.from_PLASMA(pla)
    for itn in range(100):
        eergy.solve_Te(pla, dt)
        
    eergy.to_PLASMA(pla)
    pla.update_plasma()
    
    mesh.plot_var(var=[pla.Te, pla.Ti], 
                  var_name=['E Temperature', 'Ion Temperature'],
                  fname='Prerun_Temperature.png')
    ########################################
    

    ##########################################
    ########## main loop for plasma ##########
    ##########################################
    ne_ave, ni_ave, Te_ave = [], [], []  
    time = []
   
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
        