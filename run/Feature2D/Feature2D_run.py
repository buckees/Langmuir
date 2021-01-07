"""Feature Model 2D. Main program."""

from packages.Model.Common.Yaml import PARAMETER
from packages.Model.Common.Particle import PARTICLE
from packages.Model.Feature2D.Feature2D_mesh import MESH2D
from packages.Model.Feature2D.Feature2D_rflct import REFLECT
from packages.Model.Feature2D.Feature2D_rct import REACT
from packages.Model.Feature2D.Feature2D_main import MAIN

# init operation parameters
oper = PARAMETER()
oper.num_ptcl = 10000
oper.max_step = 1000
oper.step_fac = 0.5
oper.max_rflct = 5
oper.surf_norm_range = 3
oper.surf_norm_mode = 'Sum Vector'
oper.num_plot = 5
oper.prob_rflct = 0.5

# init mesh
mesh = MESH2D()
# readin mesh
fname = 'SiEtch_Base_Mesh'
mesh.readin_mesh(fname)

# init ptcl
ptcl = PARTICLE()
ptcl.customize_ptcl('Ion', 40, 1)

# init reflection
rflct = REFLECT()

# init reaction
rct = REACT()

MAIN(oper, ptcl, mesh, rct, rflct)
