"""Feature Model 2D. Main program."""

from packages.Model.Common.Particle import PARTICLE
from packages.Model.Sheath2D.Sheath2D_ops import PARAMETER
from packages.Model.Sheath2D.Sheath2D_main import MAIN

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 10000
oper.max_step = 1000
oper.d_sh = 0.01
oper.wfr_loc = 0.0
oper.dt = 1e-8
oper.Vdc = 100.0
oper.Vrf = 50.0
oper.freq = 1

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

MAIN(oper, ptcl)
