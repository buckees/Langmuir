"""Sheath Model 2D. Main program."""

import numpy as np
from math import sin

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Sheath2D.Sheath2D_main import MAIN
from packages.Model.Sheath2D.Sheath2D_field import FIELD_SHEATH
from packages.Model.Common.Particle_Mover import EULER_MOVE, LEAPFROG
from packages.Constants import PI

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 10000
oper.max_step = 1000
oper.d_sh = 0.01
oper.wfr_loc = 0.0
oper.dt = 1e-8
oper.Vdc = 100.0
oper.Vrf = 20.0
oper.freq = 1

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

field = FIELD_SHEATH('Sheath')
field.add_para(d_sh=oper.d_sh, Vdc=oper.Vdc, Vrf=oper.Vrf, 
               freq=oper.freq, phi=0.0)

MAIN(oper, ptcl, field, move=EULER_MOVE)
