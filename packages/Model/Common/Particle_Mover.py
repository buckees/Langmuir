"""
Particle_Mover.py solves the Newton's equation of motion,
    dx/dt = v; dv/dt = F/m
    takes obj as input - PARTICLE(), MULTI_PARTICLE(), FIELD()
"""

import numpy as np

from packages.Constants import (AMU, UNIT_CHARGE)

def EULER_MOVE(ptcl, field, dt):
    """
    Solve the Newton's equatoin of motion by Euler method.
    
    x(t + dt) = x(t) + v(t) * dt
    v(t + dt) = v(t) + F(t)/m * dt
    x: arr(3) of float, position in (x, z, y)
    v: arr(3) of float, velocity in (x, z, y)
    ptcl: PARTICLE() or MULTI_PARTICLE() obj
    field: FIELD() obj
    """
    posn = ptcl.posn + ptcl.vel*dt
    vel = ptcl.vel + field.E*(ptcl.charge*UNIT_CHARGE)/(ptcl.mass*AMU)*dt
    return posn, vel

def S_EULER(ptcl, field, dt):
    """
    Semi-explicit Euler: 1st order symplectic integrator.
    
    v(t + dt) = v(t) + F(t)/m * dt
    x(t + dt) = x(t) + v(t + dt) * dt
    """    
    vel = ptcl.vel + field.E*(ptcl.charge*UNIT_CHARGE)/(ptcl.mass*AMU)*dt
    posn = ptcl.posn + ptcl.vel*dt
    return posn, vel