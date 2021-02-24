"""Examples."""

import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

from packages.Mesh.Geom import (FEATMOD2D, DOMAIN2D, 
                                RECTANGLE, TRIANGLE, CIRCLE)
from packages.Mesh.Mesh import MESH2D

# build the geometry
Feat2d = FEATMOD2D(name='Si_Etch_v02', is_cyl=False)
#               (left, bottom), (width, height)
domain2d = DOMAIN2D((0.0, 0.0),    (200.0e-9, 500.0e-9))
Feat2d.add_domain(domain2d)

# Add metal wall to all boundaries
# In Metal, vector potential A = 0
#                        (left, bottom), (right, top)
SiO2 = RECTANGLE('SiO2_', (0.0e-9, 0.0e-9), (200.0e-9, 50.0e-9))
Feat2d.add_shape(SiO2)
Si = RECTANGLE('Si_', (0.0e-9, 50.0e-9), (200.0e-9, 400.0e-9))
Feat2d.add_shape(Si)
PR = RECTANGLE('PR_', (0.0e-9, 400.0e-9), (200.0e-9, 450.0e-9)) 
Feat2d.add_shape(PR)

# Dig a hole
Vac_Rect = RECTANGLE('Plasma', (75.0e-9, 200.0e-9), (125.0e-9, 500.0e-9)) 
Feat2d.add_shape(Vac_Rect)
Vac_Circ = CIRCLE('Plasma', (100.0e-9, 200.0e-9), 60.0e-9)
Feat2d.add_shape(Vac_Circ)

# Cut a sliding mask
Vac_Trgl_1 = TRIANGLE('Plasma', (50.0e-9, 450.0e-9), (75.0e-9, 450.0e-9), 
                    (75.0e-9, 400.0e-9))
Feat2d.add_shape(Vac_Trgl_1)
Vac_Trgl_2 = TRIANGLE('Plasma', (125.0e-9, 450.0e-9), (140.0e-9, 450.0e-9), 
                    (125.0e-9, 400.0e-9))
Feat2d.add_shape(Vac_Trgl_2)
Vac_Trgl_3 = TRIANGLE('Plasma', (100.0e-9, 400.0e-9), (100.0e-9, 200.0e-9), 
                    (25.0e-9, 200.0e-9))
Feat2d.add_shape(Vac_Trgl_3)


Feat2d.plot(figsize=(4, 4), ihoriz=1)
print(Feat2d)

# generate mesh to imported geometry
mesh2d = MESH2D(import_geom=Feat2d)
mesh2d.gen_mesh(ngrid=(50, 125))
mesh2d.plot(figsize=(4, 4))
