#!/usr/bin/env python3

# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 20
#
#  STEP import and manipulation, geometry partitioning
#
# ------------------------------------------------------------------------------

# The OpenCASCADE CAD kernel allows to import STEP files and to modify them. In
# this tutorial we will load a STEP geometry and partition it into slices.

import gmsh
import math
import os
import sys

gmsh.initialize()

gmsh.model.add("chip")

# Load a STEP file (using `importShapes' instead of `merge' allows to directly
# retrieve the tags of the highest dimensional imported entities):
chipSolid = gmsh.model.occ.importShapes('Tilting-Chip-simplified.stp')

# Get the bounding box of the volume:
xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(
    chipSolid[0][0], chipSolid[0][1]
)
# Create a box and cut out the geometry (to obtain only the void space)
boundingBox = gmsh.model.occ.addBox(xmin, ymin, zmin, xmax-xmin, ymax-ymin, zmax-zmin)
chipVoid = gmsh.model.occ.cut([(3, boundingBox)], chipSolid)
gmsh.model.occ.synchronize()

entities = gmsh.model.getEntities()
for e in entities:
    print(e)

gmsh.write("chip.step")

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
