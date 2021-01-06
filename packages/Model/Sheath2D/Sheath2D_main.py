"""Sheath Model 2D. Main program."""

import numpy as np
from math import sin
import matplotlib.pyplot as plt

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Common.Field import FIELD
from packages.Model.Common.Particle_Mover import EULER_MOVE
from packages.Constants import (PI, AMU, UNIT_CHARGE)

# Physics parameters
w_loc = 0.0 # by default

dsh = 0.01
dt = 1e-8
Vdc, Vrf = 100, 20
# freq = 13.56e6 # T = 70 ns
freq = 1

def Efunc(Vdc, Vrf, dsh, freq, phi0, t):
    """Calc sheath E-field."""
    E = np.zeros(3)
    E[1] = -Vdc/dsh - Vrf/dsh*sin(2*PI*freq*t + phi0)
    return E

# Model parameters
num_ptcl = 10000
max_step = 10000

ptcl = PARTICLE('M+')
ptcl.customize_ptcl('Ion', 40, 1)

field = FIELD('Sheath')
field.add_Efunc(Efunc)



def MAIN(oper, ptcl, coll=None):
    """
    MAIN() actually runs the feature model.
    oper: OPERATION(obj), contains all operation parameters.
    ptcl: PARTICLE(obj), contains all particle information.
    coll: COLLISION(obj), contains all collision information.
    """

    erg = list()
    for i in range(oper.num_ptcl):
        
        ########## init ptcl ##########
        ptcl.update_state(True)
        ptcl.update_posn(np.array([0.0, oper.d_sh, 0.0]))
        ptcl.update_vel(np.zeros(3))
        ################################
        
        dt = oper.dt
        t = 0
        step = 0
        phi0 = np.random.uniform(0.0, 2*PI)
        
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
            field.update_E(field.Efunc(oper.Vdc, oper.Vrf, oper.d_sh, 
                                       oper.freq, phi0, t))
            posn_next, vel_next = EULER_MOVE(ptcl, field, dt)
            ptcl.update_posn(posn_next)
            ptcl.update_vel(vel_next)
            t += dt
            step += 1
            if ptcl.posn[1] < oper.wfr_loc:
                ptcl.update_state(False)
                erg.append(ptcl.vel2erg())
            if step > max_step:
                ptcl.update_state(False)
    
    ########## plot results ##########
    print(f'{num_ptcl} particles are launched.' 
          + f'\n{len(erg)} particles are collected by the wafer.')
    
    
    fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=600,
                           constrained_layout=True)
    ax.hist(erg, 2*int(Vrf), density=True)
    ax.set_xlim(40, 160)
    ax.set_title('Ion Energy Distribution')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Count')
    fig.savefig(f'Vdc{oper.Vdc}_Vrf{oper.Vrf}.png', dpi=600)
    plt.close()
