"""Feature Model 2D. Main program."""

import numpy as np
from math import sin, cos
from copy import deepcopy
import random
import pandas as pd

def MAIN(oper, ptcl, mesh, chem, rct, rflct):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    ptcl: PARTICLE(obj), contains all particle information.
    mesh: MESH2D(obj), contains all mesh informatin.
    chem: CHEMISTRY(obj), contains all reaction information.
    rct: REACTION(obj), contains all reaction informatino.
    rflct: REFLECTION(obj), contains all reflection information.
    """
        
    # init diagnostics
    if oper.idiag:
        init_posn, init_vel = list(), list()
        rec_traj, rec_surf, rec_mesh = [], [], []
    
    # update mesh dict
    num_mat = len(mesh.dict_num2mat)
    for mat in chem.df_mat['Material']:
        if mat not in mesh.dict_mat2num.keys():
            num_mat += 1
            mesh.dict_mat2num[mat] = num_mat
            mesh.dict_num2mat[num_mat] = mat
    # readin flux
    df_flux = pd.read_csv(oper.fname + '_Flux.csv', header=0)
    sp_list = df_flux['Species'].tolist()
    weight = df_flux['Flux'].tolist()
    
    delta_L = (mesh.res*oper.step_fac).min()

    for k in range(oper.num_ptcl):
        
        ########## choose particle #########
        chosen_ptcl = random.choices(sp_list, weights=weight, k=1)[0]
        ptcl.select_ptcl(chosen_ptcl)
        ####################################
        
        idiag = False
        if k > oper.num_ptcl - 20:
            idiag = True
        
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
            init_posn.append(ptcl.posn)
            init_vel.append(ptcl.vel)
        ################################
        
        ########## record initial position ##########
        if idiag:
            rec_traj.append([])
            rec_traj[-1].append(ptcl.posn.copy())
        #############################################
        
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
            if not (mesh.left < ptcl.posn[0] < mesh.right):
                posn = deepcopy(ptcl.posn)
                posn[0]= mesh.left + ((posn[0] - mesh.left) % mesh.width)
                ptcl.update_posn(posn)
            ################################
            
            # check if the ptcl is dead
            if not ptcl.isAlive:
                # record ptcl posn when dead
                if idiag:
                    rec_traj[-1].append(deepcopy(ptcl.posn))
                break
            
            # check hit
            hit_mat, hit_idx = mesh.check_hit(ptcl.posn[0:2])
            if hit_mat:
                # record the hit point
                if idiag:
                    rec_traj[-1].append(deepcopy(ptcl.posn))
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
                        break
                    
                    vel = deepcopy(ptcl.vel)
                    vel[0:2] = rflct.rflct(vel[0:2])
                    ptcl.update_vel(vel)
                    num_rflct += 1
                else:
                    # now ireact = 1
                    mesh.change_mat(hit_idx)
                    ptcl.update_state(False)
            # check if the ptcl is dead
            if not ptcl.isAlive:
                # record ptcl posn when dead
                if idiag:
                    rec_traj[-1].append(deepcopy(ptcl.posn))
                break
    
        if idiag:
            rec_traj[-1] = np.array(rec_traj[-1]).T
    
    np.save('Feat2D_initPosn.npy', init_posn)
    np.save('Feat2D_initVel.npy', init_vel)
   
    # rec_surf = []
    # for temp_idx in mesh.surf_set:
    #     temp_svec, temp_stheta = mesh.calc_surf_norm(temp_idx, 
    #                                                   radius=surf_norm_range, 
    #                                                   imode=surf_norm_mode)
    #     rec_surf.append([temp_idx, temp_svec])
    
    # colMap = copy(cm.Accent)
    # colMap.set_under(color='white')
    
    # def plot_traj(ax, traj):
    #     # for i in range(10):
    #     #     ax.plot(traj[num_ptcl - i - 1][0, :], traj[num_ptcl - i - 1][1, :],
    #     #             marker='o', markersize=0.3, linestyle='-', linewidth=0.1)
    #     for temp_ptcl in traj:
    #         ax.plot(temp_ptcl[0, :], temp_ptcl[1, :],
    #                 marker='o', markersize=0.3, linestyle='-', linewidth=0.1)
    
    # def plot_surf_norm(ax, posn, svec):
    #     ax.quiver(posn[0], posn[1],
    #               svec[0], svec[1], 
    #               # scale=50, units='xy',
    #               # headwidth=1, headlength=1, lw=0.01, edgecolors='k',
    #               width=0.001)
    
    
    # def plot_mesh(mat, ith):
    #     fig, axes = plt.subplots(1, 2, figsize=(16, 8),
    #                               constrained_layout=True)
        
    #     ax = axes[0]
    #     ax.contourf(mesh.x, mesh.z, mat, cmap=colMap, vmin=0.2, extend='both')
    #     ax.set_xlim(0.0, mesh.width)
    #     ax.set_ylim(0.0, mesh.height)
    #     plot_traj(ax, rec_traj)
        
    #     ax = axes[1]
    #     ax.scatter(mesh.x, mesh.z, c=mat, s=1, cmap=colMap, vmin=0.2)
    #     ax.set_xlim(0.0, mesh.width)
    #     ax.set_ylim(0.0, mesh.height)
    #     plot_traj(ax, rec_traj)
    #     for item in rec_surf:
    #         temp_idx, temp_svec = item
    #         temp_posn = np.array([mesh.x[temp_idx], mesh.z[temp_idx]])
    #         plot_surf_norm(ax, temp_posn, temp_svec)
        
    #     fig.savefig('mat_%d.png' % ith, dpi=300)
    
    # for ith, mat in enumerate(rec_mesh):
    #     plot_mesh(mat, ith+1)
