"""Sheath Model 2D. Main program."""

import numpy as np
from math import exp
from copy import deepcopy
import pandas as pd

from packages.Constants import PI

def MAIN(oper, ptcl, field, coll, move, stats=None):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    ptcl: PARTICLE(obj), contains all particle information.
    field: FIELD(obj), contains all field information.
    coll: COLLISION(obj), contains all collision information.
    move: Func, selected from Particle_Mover
    stats: STATS(obj), contains all diagnostic information.
    """

    vel = list()
    for i in range(oper.num_ptcl):
        ########## init ptcl ##########
        ptcl.update_state(True)
        ptcl.init_posn()
        ptcl.init_vel_vFunc()
        if oper.idiag:
            # stats.append_row()
            # print(stats.df)
            # print(stats.df.iloc[-1])
            # print(stats.df.iloc[-1].loc['Collision'])
            stats.df.loc[i,'Collision'] = 0
            ptcl_erg, ptcl_ang = ptcl.vel2erg()
            stats.df.loc[i, ['Init_Vel_x', 'Init_Vel_z', 'Init_Vel_y']] = \
                ptcl.vel
            stats.df.loc[i, 'Init_Erg'] = deepcopy(ptcl_erg)
            stats.df.loc[i, 'Init_Ang'] = deepcopy(ptcl_ang)
        ################################
        
        dt = oper.dt
        step = 0
        phi0 = np.random.uniform(0.0, 2*PI)
        t = phi0/(2*PI*oper.freq)
        
        ###############################################
        ########## main loop for pctl launch ##########
        ###############################################
        while ptcl.isAlive:
            # make sure not over shoot
            if ptcl.vel[1]:
                dt1 = (oper.wfr_loc - ptcl.posn[1])/ptcl.vel[1]
                if dt1 < dt:
                    dt = dt1*1.001
            # move
            field.update_E(field.calc_E(t))
            posn_next, vel_next = move(ptcl, field, t, dt)
            ptcl.update_posn(posn_next)
            ptcl.update_vel(vel_next)
            t += dt
            step += 1
            
            if ptcl.posn[1] < oper.wfr_loc:
                ptcl.update_state(False)
                vel.append(deepcopy(ptcl.vel))
                if oper.idiag:
                    stats.df.loc[i, 'hitWafer'] = True
            if step > oper.max_step:
                ptcl.update_state(False)
                if oper.idiag:
                    stats.df.loc[i, 'hitWafer'] = False

            # collision
            coll_freq = coll.func_CollFreq(ptcl.vel)
            prob_coll = 1.0 - exp( - coll_freq * dt)
            rand = np.random.uniform(0.0, 1.0)
            if rand < prob_coll:
                vel_new = coll.func_ReinitVel(ptcl.vel)
                ptcl.update_vel(vel_new)
                if oper.idiag:
                    stats.df.loc[i, 'Collision'] += 1
        
        # now particle is marked as dead
        if oper.idiag:
            ptcl_erg, ptcl_ang = ptcl.vel2erg()
            stats.df.loc[i, ['End_Vel_x', 'End_Vel_z', 'End_Vel_y']] = \
                ptcl.vel
            stats.df.loc[i, 'End_Erg'] = deepcopy(ptcl_erg)
            stats.df.loc[i, 'End_Ang'] = deepcopy(ptcl_ang)
    ########## plot results ##########
    print(f'{oper.num_ptcl} particles are launched.' 
          + f'\n{len(vel)} particles are collected by the wafer.')
    
    if oper.idiag:
        return vel, stats
    else:
        return vel