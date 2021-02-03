"""Test Geom2d."""

import os
import glob
for i in glob.glob("*.png"):
    os.remove(i)

from packages.Mesh.Geom import RCTMOD2D, DOMAIN2D, RECTANGLE
from packages.Mesh.Mesh import MESH2D

# build the geometry
ICP2d = RCTMOD2D(name='ICP2D', is_cyl=False)
#               (left, bottom), (width, height)
domain2d = DOMAIN2D((-0.25, 0.0),    (0.5, 0.4))
ICP2d.add_domain(domain2d)

# Add metal wall to all boundaries
# In Metal, vector potential A = 0
#                        (left, bottom), (right, top)
top = RECTANGLE('Metal', (-0.25, 0.38), (0.25, 0.4))
ICP2d.add_shape(top)
bott = RECTANGLE('Metal', (-0.25, 0.0), (0.25, 0.02))
ICP2d.add_shape(bott)
left = RECTANGLE('Metal', (-0.25, 0.0), (-0.23, 0.4)) 
ICP2d.add_shape(left)
right = RECTANGLE('Metal', (0.23, 0.0), (0.25, 4.0))
ICP2d.add_shape(right)
ped = RECTANGLE('Metal', (-0.20, 0.0), (0.20, 0.1))
ICP2d.add_shape(ped)


# Add quartz to separate coil area and plasma area
# Quartz conductivity = 1e-5 S/m
quartz = RECTANGLE('Quartz', (-0.23, 0.3), (0.23, 0.32))
ICP2d.add_shape(quartz)

# Add air to occupy the top coil area to make it non-plasma
# Air concudctivity = 0.0 S/m
air = RECTANGLE('Air', (-0.23, 0.32), (0.23, 0.38))
ICP2d.add_shape(air)

# Add coil within air and overwirte air
# coil 1, 2, 3: J = -J0*exp(iwt)
# coil 4, 5, 6: J = +J0*exp(iwt)
coil1 = RECTANGLE('Coil', (-0.20, 0.34), (-0.18, 0.36))
ICP2d.add_shape(coil1)
coil2 = RECTANGLE('Coil', (-0.14, 0.34), (-0.12, 0.36))
ICP2d.add_shape(coil2)
coil3 = RECTANGLE('Coil', (-0.08, 0.34), (-0.06, 0.36))
ICP2d.add_shape(coil3)
coil4 = RECTANGLE('Coil', (0.18, 0.34), (0.20, 0.36))
ICP2d.add_shape(coil4)
coil5 = RECTANGLE('Coil', (0.12, 0.34), (0.14, 0.36))
ICP2d.add_shape(coil5)
coil6 = RECTANGLE('Coil', (0.06, 0.34), (0.08, 0.36))
ICP2d.add_shape(coil6)
 
ICP2d.plot(figsize=(10, 4), ihoriz=1)
print(ICP2d)

# generate mesh to imported geometry
mesh2d = MESH2D(import_geom=ICP2d)
mesh2d.gen_mesh(ngrid=(50, 40))
mesh2d.plot(figsize=(10, 4))
