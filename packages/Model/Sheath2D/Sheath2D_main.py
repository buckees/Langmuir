"""Sheath Model 2D. Main program."""

import numpy as np
import matplotlib.pyplot as plt

from packages.Constants import PI

def plot_erg(erg):
    """Plot energy distribution."""
    fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                           constrained_layout=True)
    ax.hist(erg, bins=50, density=True)
    ax.set_title('Ion Energy Distribution')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Count')
    fig.savefig('erg_distrb.png', dpi=600)
    plt.close()


def MAIN(oper, ptcl, field, coll=None, move=None):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    ptcl: PARTICLE(obj), contains all particle information.
    field: FIELD(obj), contains all field information.
    coll: COLLISION(obj), contains all collision information.
    move: Func, selected from Particle_Mover
    """

    erg, ang = list(), list()
    for i in range(oper.num_ptcl):
        
        ########## init ptcl ##########
        ptcl.update_state(True)
        ptcl.update_posn(np.array([0.0, oper.d_sh, 0.0]))
        ptcl.update_vel(np.zeros(3))
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
                ptcl_erg, ptcl_ang = ptcl.vel2erg()
                erg.append(ptcl_erg)
                ang.append(ptcl_ang)
            if step > oper.max_step:
                ptcl.update_state(False)
    
    ########## plot results ##########
    print(f'{oper.num_ptcl} particles are launched.' 
          + f'\n{len(erg)} particles are collected by the wafer.')
    
    if oper.iplot:
        plot_erg(erg)
    
    return erg, ang