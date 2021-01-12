"""
Particle_Mover.py solves the Newton's equation of motion,
    dx/dt = v; dv/dt = F/m
    takes obj as input - PARTICLE(), MULTI_PARTICLE(), FIELD()
"""

import numpy as np

from packages.Constants import (AMU, UNIT_CHARGE)

def EULER_MOVE(ptcl, field, t, dt):
    """
    Solve the Newton's equatoin of motion by Euler method.
    
    x(t + dt) = x(t) + v(t) * dt
    v(t + dt) = v(t) + F(t)/m * dt
    posn: arr(3) of float, position in (x, z, y)
    vel: arr(3) of float, velocity in (x, z, y)
    ptcl: PARTICLE() or MULTI_PARTICLE() obj
    field: FIELD() obj
    t: unused
    dt: float, unit in s, timestep
    """
    Coulomb_const = (ptcl.charge*UNIT_CHARGE)/(ptcl.mass*AMU)
    posn = ptcl.posn + ptcl.vel*dt
    vel = ptcl.vel + field.calc_E(t)*Coulomb_const*dt
    return posn, vel

def S_EULER(ptcl, field, t, dt):
    """
    Semi-explicit Euler: 1st order symplectic integrator.
    
    v(t + dt) = v(t) + F(t)/m * dt
    x(t + dt) = x(t) + v(t + dt) * dt
    posn: arr(3) of float, position in (x, z, y)
    vel: arr(3) of float, velocity in (x, z, y)
    ptcl: PARTICLE() or MULTI_PARTICLE() obj
    field: FIELD() obj
    t: unused
    dt: float, unit in s, timestep
    """
    Coulomb_const = (ptcl.charge*UNIT_CHARGE)/(ptcl.mass*AMU)
    vel = ptcl.vel + field.calc_E(t)*Coulomb_const*dt
    posn = ptcl.posn + ptcl.vel*dt
    return posn, vel

def LEAPFROG(ptcl, field, t, dt):
    """
    Semi-explicit Euler: 1st order symplectic integrator.
    
    v(t + dt) = v(t) + F(t)/m * dt
    x(t + dt) = x(t) + v(t + dt) * dt
    posn: arr(3) of float, position in (x, z, y)
    vel: arr(3) of float, velocity in (x, z, y)
    ptcl: PARTICLE() or MULTI_PARTICLE() obj
    field: FIELD() obj
    t: float, unit in s, current time
    dt: float, unit in s, timestep
    """
    Coulomb_const = (ptcl.charge*UNIT_CHARGE)/(ptcl.mass*AMU)
    Acur = field.calc_E(t)*Coulomb_const
    # CWong's note: Efunc interface would be better if only need to poass time variable
    Anew = field.calc_E(t + dt) * Coulomb_const
    posn = ptcl.posn + ptcl.vel*dt + 0.5*Acur*dt*dt
    vel = ptcl.vel + 0.5*(Acur + Anew)*dt
    return posn, vel
    