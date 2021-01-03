"""Feature Model 2D. Main program."""

import numpy as np
from math import sin, cos
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from copy import copy, deepcopy


from packages.Model.Feature2D.Feature2D_ops import (width, height, 
                                                res_x, res_z, num_ptcl, ibc, 
                          threshold, max_rflct, idstrb, step_fac, max_step,
                          num_plot, surf_norm_range, surf_norm_mode)
from packages.Model.Common.Particle import PARTICLE
from packages.Model.Feature2D.Feature2D_mesh import MESH2D
from packages.Model.Feature2D.Feature2D_rflct import REFLECT

# init MESHGRID obj
mesh = MESH2D()
# readin mesh
fname = 'SiEtch_Base_Mesh'
mesh.readin_mesh(fname)

delta_L = (mesh.res*step_fac).min()

# species information is imported from species
# Initialize the PARTICLE() object
ptcl = PARTICLE('Ar+')
ptcl.customize_ptcl('Ion', 40, 1)
ptcl_rflct = REFLECT()

# init diagnostics
rec_traj, rec_surf, rec_mesh = [], [], []

num_ptcl = 10000
for k in range(num_ptcl):
    if (k + 1) % int(num_ptcl/num_plot) == 0:
        print('%d particles are launched!' % (k+1))
        mesh.plot(dpi=300, fname='nptcl=%d.png' % (k+1))
        # rec_mesh.append(deepcopy(mesh.mat))
    ptcl.isAlive = True
    # gen rand position at the top, in (x, z, y)
    posn = mesh.init_ptcl_posn()
    # pass posn to ptcl
    ptcl.update_posn(posn)
    # record initial position
    if k > num_ptcl - 20:
        rec_traj.append([])
        rec_traj[-1].append(ptcl.posn.copy())
    
    # gen rand velocity
    vel = np.zeros(3)
    mu, sigma = 0.0, 0.1  # default mean and standard deviation
    theta = np.random.normal(mu, sigma)
    vel[0], vel[1] = sin(theta), -cos(theta)
    # pass vel to ptcl
    ptcl.update_vel(vel)

    num_rflct = 0
#    while imove_ptcl == 1 and num_rflct < 5:
    for i in range(max_step):
        # move the ptcl in space by delta_L
        ptcl.move_in_space(delta_L)
        # periodic b.c. at left and right bdry, make it a func?
        if ptcl.posn[1] >= (mesh.top - mesh.dz*0.5):
            ptcl.update_state(False)
        if not (mesh.left < ptcl.posn[0] < mesh.right):
            posn = deepcopy(ptcl.posn)
            posn[0]= mesh.left + ((posn[0] - mesh.left) % mesh.width)
            ptcl.update_posn(posn)
        
        # check if the ptcl is dead
        if not ptcl.isAlive:
            # record ptcl posn when dead
            if k > num_ptcl - 20:
                rec_traj[-1].append(deepcopy(ptcl.posn))
            break
        
        hit_mat, hit_idx = mesh.check_hit(ptcl.posn[0:2])
        if hit_mat:
            # record the hit point
            if k > num_ptcl - 20:
                rec_traj[-1].append(deepcopy(ptcl.posn))
            # at this position, th ptcl hits a mat
            mat_name = mesh.mat_dict[hit_mat]
            # calc surf norm
            ptcl_rflct.svec, ptcl_rflct.stheta = \
                mesh.calc_surf_norm(hit_idx, radius=surf_norm_range, 
                                    imode=surf_norm_mode)
            # decide wehter a reflection or reaction
            rand = np.random.uniform(0.0, 1.0)
            rflct = REFLECT(ptcl.name, mat_name, 1.0)
            prob = rflct.calc_prob()
            if rand < prob:
                # check max rflct
                if num_rflct > max_rflct:
                    ptcl.update_state(False)
                    break
                
                vel = deepcopy(ptcl.vel)
                temp_vel = ptcl_rflct.rflct(vel[0:2])
                vel[0] = temp_vel[0]
                vel[1] = temp_vel[1]
                ptcl.update_vel(vel)
                num_rflct += 1
            else:
                # now ireact = 1
                mesh.update_mat(hit_idx, threshold)
                # mesh.find_float_cell()
                ptcl.update_state(False)
        # check if the ptcl is dead
        if not ptcl.isAlive:
            # record ptcl posn when dead
            if k > num_ptcl - 20:
                rec_traj[-1].append(deepcopy(ptcl.posn))
            break

    if k > num_ptcl - 20:
        rec_traj[-1] = np.array(rec_traj[-1]).T


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
