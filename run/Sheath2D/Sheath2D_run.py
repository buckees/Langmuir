"""Sheath Model 2D. Main program."""

import numpy as np
from math import sin

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Common.Field import FIELD
from packages.Model.Sheath2D.Sheath2D_ops import PARAMETER
from packages.Model.Sheath2D.Sheath2D_main import MAIN
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

# init field
def Efunc(Vdc, Vrf, dsh, freq, phi0, t):
    """Calc sheath E-field."""
    E = np.zeros(3)
    E[1] = -Vdc/dsh - Vrf/dsh*sin(2*PI*freq*t + phi0)
    return E

field = FIELD('Sheath')
field.add_Efunc(Efunc)

MAIN(oper, ptcl, field)
