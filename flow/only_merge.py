#!/usr/bin/env python3

import time
import numpy as np
from mgtools import MGread, ScalarField, VectorField, h5Viewer, vtkViewer

Nstart = 0
Nstop = 1131
    
x = np.linspace( -8.0, 40.0, 2400+1)
y = np.linspace( -8.0, 8.0, 800+1)
z = np.linspace(  4.0, 8.0, 200+1)
lambda2 = ScalarField(x=x, y=y, z=z)
velocity = VectorField(x=x, y=y, z=z)
vorticity = VectorField(x=x, y=y, z=z)
    
viewer = h5Viewer('cyl_34_1.h5')
    
for grid in range(Nstart, Nstop):
    tic = time.perf_counter()

    viewer.set_group('/grid-{}'.format(grid))
    lambda2Grid = viewer.read("lambda2")
    velocityGrid = viewer.read("velocity")
    vorticityGrid = viewer.read("vorticity")

    # Strip = 2, get rid of interface lines for grid topology improvement
    # Strip = 1, to get rid of the 'grids' on the surface of the iso-surfaces
    lambda2.merge(lambda2Grid, strip=1)
    velocity.merge(velocityGrid, strip=1)
    vorticity.merge(vorticityGrid, strip=1)

    toc = time.perf_counter()
    print("Merged grid {0} in {1:.3f} sec.".format(grid, toc-tic))
    
viewer.close()

magVorticity = vorticity.mag()
    
# Save data to disk as VTK
viewer = vtkViewer('cyl_34_1.vtk', equidistant=True)
viewer.write(lambda2, "lambda2")
viewer.write(velocity, "velocity")
viewer.write(vorticity, "vorticity")
viewer.write(magVorticity, "magVorticity")
viewer.close()
