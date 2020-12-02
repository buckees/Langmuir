"""Examples."""

import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

from packages.Mesh.LngmrMod_Geom import (FeatMod2D, Domain2D, 
                                         Rectangle, Circle)
from packages.Mesh.LngmrMod_Mesh import Mesh2D

# build the geometry
Feat2d = FeatMod2D(name='Feat2D', is_cyl=False)
#               (left, bottom), (width, height)
domain2d = Domain2D((0.0, 0.0),    (200.0e-9, 500.0e-9))
Feat2d.add_domain(domain2d)

# Add metal wall to all boundaries
# In Metal, vector potential A = 0
#                        (left, bottom), (right, top)
SiO2 = Rectangle('SiO2_', (0.0e-9, 0.0e-9), (200.0e-9, 50.0e-9))
Feat2d.add_shape(SiO2)
Si = Rectangle('Si_', (0.0e-9, 50.0e-9), (200.0e-9, 400.0e-9))
Feat2d.add_shape(Si)
PR = Rectangle('PR_', (0.0e-9, 400.0e-9), (200.0e-9, 450.0e-9)) 
Feat2d.add_shape(PR)

# Dig a hole
Vac_Rect = Rectangle('Plasma', (75.0e-9, 200.0e-9), (125.0e-9, 500.0e-9)) 
Feat2d.add_shape(Vac_Rect)
Vac_Circ = Circle('Plasma', (100.0e-9, 200.0e-9), 60.0e-9)
Feat2d.add_shape(Vac_Circ)

Feat2d.plot(figsize=(4, 4), ihoriz=1)
print(Feat2d)

# generate mesh to imported geometry
mesh2d = Mesh2D(import_geom=Feat2d)
mesh2d.gen_mesh(ngrid=(50, 125))
mesh2d.plot(figsize=(4, 4), ihoriz=1, s_size=1)
