"""Sheath Model 2D. Main program."""

import numpy as np
from math import exp
from copy import deepcopy

from packages.Constants import PI

def MAIN(oper, ptcl, field, coll, move=None):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    ptcl: PARTICLE(obj), contains all particle information.
    field: FIELD(obj), contains all field information.
    coll: COLLISION(obj), contains all collision information.
    move: Func, selected from Particle_Mover
    """

    vel = list()
    if oper.idiag:
        init_erg, init_ang = list(), list()
        erg, ang = list(), list()
    for i in range(oper.num_ptcl):
        
        ########## init ptcl ##########
        ptcl.update_state(True)
        ptcl.init_posn()
        ptcl.setVel_vFunc()
        if oper.idiag:
            ptcl_erg, ptcl_ang = ptcl.vel2erg()
            init_erg.append(ptcl_erg)
            init_ang.append(ptcl_ang)
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
                ptcl_erg, ptcl_ang = ptcl.vel2erg()
                erg.append(ptcl_erg)
                ang.append(ptcl_ang)
            if step > oper.max_step:
                ptcl.update_state(False)

            # collision
            coll_freq = coll.func_CollFreq(ptcl.vel)
            prob_coll = 1.0 - exp( - coll_freq * dt)
            rand = np.random.uniform(0.0, 1.0)
            if rand < prob_coll:
                ptcl.setVel_norm()
    
    ########## plot results ##########
    print(f'{oper.num_ptcl} particles are launched.' 
          + f'\n{len(erg)} particles are collected by the wafer.')
    
    if oper.idiag:
        return vel, erg, ang, init_erg, init_ang
    else:
        return vel