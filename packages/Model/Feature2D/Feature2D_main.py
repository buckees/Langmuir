"""Feature Model 2D. Main program."""

import numpy as np
from math import sin, cos
from copy import deepcopy
import random
import pandas as pd

def MAIN(oper, ptcl, mesh, chem, rct, rflct, stats=None):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    ptcl: PARTICLE(obj), contains all particle information.
    mesh: MESH2D(obj), contains all mesh informatin.
    chem: CHEMISTRY(obj), contains all reaction information.
    rct: REACTION(obj), contains all reaction informatino.
    rflct: REFLECTION(obj), contains all reflection information.
    stats: STATS(obj), contains all diagnostic information.
    """
            
    # update mesh dict
    num_mat = len(mesh.dict_num2mat)
    for mat in chem.df_mat['Material']:
        if mat not in mesh.dict_mat2num.keys():
            num_mat += 1
            mesh.dict_mat2num[mat] = num_mat
            mesh.dict_num2mat[num_mat] = mat
    # readin flux
    df_flux = pd.read_csv(oper.fcase + '_Flux.csv', header=0)
    sp_list = df_flux['Species'].tolist()
    weight = df_flux['Flux'].tolist()
    
    delta_L = (mesh.res*oper.step_fac).min()

    for k in range(oper.num_ptcl):
        
        ########## choose particle #########
        chosen_ptcl = random.choices(sp_list, weights=weight, k=1)[0]
        ptcl.select_ptcl(chosen_ptcl)
        if oper.idiag:
            stats.df_sp.loc[ptcl.name, 'Launched'] += 1
        ####################################
                
        ########## print progress ##########
        if (k + 1) % int(oper.num_ptcl/oper.num_plot) == 0:
            print('%d particles are launched!' % (k+1))
            mesh.plot(dpi=300, fname='nptcl=%d.png' % (k+1))
            # rec_mesh.append(deepcopy(mesh.mat))
        ####################################
        
        ########## init ptcl ##########
        # make the ptcl alive
        ptcl.update_state(True)
        ptcl.init_posn()
        if ptcl.ptype == 'Ion':
            ptcl.init_vel_vFunc()
        elif ptcl.ptype == 'Neut':
            ptcl.init_vel_norm(oper.Tn)
        else:
            print('"f{ptcl.ptype}" is not found in the database.')
        if oper.idiag:
            stats.df_sp.loc[ptcl.name, 'Init Posn'].append(deepcopy(ptcl.posn))
            stats.df_sp.loc[ptcl.name, 'Init Vel'].append(deepcopy(ptcl.vel))
            init_erg, init_ang = ptcl.vel2erg()
            stats.df_sp.loc[ptcl.name, 'Init Erg'].append(deepcopy(init_erg))
            stats.df_sp.loc[ptcl.name, 'Init Ang'].append(deepcopy(init_ang))
        ################################
        
        ###############################################
        ########## main loop for pctl launch ##########
        ###############################################
        num_rflct = 0
        for i in range(oper.max_step):
            # move the ptcl in space by delta_L
            ptcl.move_in_space(delta_L)
            
            ########## check b.c. ##########
            if ptcl.posn[1] > mesh.top:
                ptcl.update_state(False)
                if oper.idiag:
                    stats.df_sp.loc[ptcl.name, 'Escaped'] += 1
                break
            if not (mesh.left < ptcl.posn[0] < mesh.right):
                posn = deepcopy(ptcl.posn)
                posn[0]= mesh.left + ((posn[0] - mesh.left) % mesh.width)
                ptcl.update_posn(posn)
            ################################
            
            # check hit
            hit_mat, hit_idx = mesh.check_hit(ptcl.posn[0:2])
            if hit_mat:
                # at this position, th ptcl hits a mat
                # calc surf norm
                rflct.svec, rflct.stheta = \
                    mesh.calc_surf_norm(hit_idx, radius=oper.surf_norm_range, 
                                        imode=oper.surf_norm_mode)
                # determine the reaction/reflection
                ptcl_erg, ptcl_ang = ptcl.vel2erg()
                hit_mat_name = mesh.dict_num2mat[hit_mat]
                chosen_rct = chem.choose_react(ptcl.name, hit_mat_name, 
                                               ptcl_erg, 0.0)
                
                if chosen_rct['Reaction_Type'] == 'Reflect':
                    # check max rflct
                    if num_rflct > oper.max_rflct:
                        ptcl.update_state(False)
                        if oper.idiag:
                            stats.df_sp.loc[ptcl.name, 'Terminated'] += 1
                        break
                    
                    vel = deepcopy(ptcl.vel)
                    vel[0:2] = rflct.rflct(vel[0:2])
                    ptcl.update_vel(vel)
                    num_rflct += 1
                else:
                    # now ireact = 1
                    mesh.change_mat(hit_idx)
                    ptcl.update_state(False)
                    if oper.idiag:
                        stats.df_sp.loc[ptcl.name, 'Etch'] += 1
                        stats.df_mat.loc[hit_mat_name, 'Etch'] += 1
                    break
        # when ptcl is dead        
        if oper.idiag:
            stats.df_sp.loc[ptcl.name, 'Reflection'].append(deepcopy(num_rflct))
    ################ OUTPUT DIAG INFO #####################
    if oper.idiag:
        # for sp, erg in stats.erg.items():
        #     np.save('initErg_' + sp, erg)
        # for sp, ang in stats.ang.items():
        #     np.save('initAng_' + sp, ang)
        stats.df_sp['Escaped Pct'] = \
            stats.df_sp['Escaped']/stats.df_sp['Launched']
        stats.df_sp['Etch Pct'] = stats.df_sp['Etch']/stats.df_sp['Launched']
        stats.df_sp['Terminated Pct'] = \
            stats.df_sp['Terminated']/stats.df_sp['Launched']

        return stats