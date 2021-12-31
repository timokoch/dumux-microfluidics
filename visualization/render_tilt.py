#!/usr/bin/env pvpython

from paraview.simple import *
import numpy as np
import time
import glob

# create a new 'Legacy VTK Reader'
filename = "synth.vtk" # "chip_void.vtk"
chip_voidvtk = LegacyVTKReader(
    registrationName='chip_void',
    FileNames=['/Users/pumbaa/Sync/Data/LabOnAChip/' + filename]
)

files = [f"/Users/pumbaa/dune-master/dumux-microfluidic/build-cmake/test/reservoir/intersections-{i}.vtu" for i in range(36)]
intersections = XMLUnstructuredGridReader(registrationName='intersections', FileName=files)

animationScene = GetAnimationScene()
animationScene.UpdateAnimationUsingDataTimeSteps()
intersections.TimeArray = 'None'

# get active view
renderView = GetActiveViewOrCreate('RenderView')

# translate to origin (some how the geometry is not centered)
transform = Transform(Input=chip_voidvtk)
if "synth" in filename:
    transform.Transform.Translate = [0.0, 0.0, 0.0]
else:
    transform.Transform.Translate = [7.3, -2.0, 3.75]

transform2 = Transform(Input=transform)
transform2.Transform.Rotate = [0.0, 0.0, 0.0]

Show()
disp = GetDisplayProperties(transform2, view=renderView)
disp.Opacity = 0.2

transformWater = Transform(Input=intersections)
transformWater.Transform.Rotate = [0.0, 0.0, 0.0]

Show()
dispWater = GetDisplayProperties(transformWater, view=renderView)
dispWater.AmbientColor = [0.6666666666666666, 1.0, 1.0]
dispWater.DiffuseColor = [0.6666666666666666, 1.0, 1.0]

Render()

renderView = GetActiveViewOrCreate('RenderView')

cylinder1 = Cylinder(registrationName='Cylinder1')

# Properties modified on cylinder1
cylinderHalfHeight = 5.0
cylinder1.Resolution = 30
cylinder1.Height = cylinderHalfHeight*2
cylinder1.Radius = 0.5

# show data in view
cylinder1Display = Show()
renderView.Update()

cylinder1Display.Orientation = [90.0, 0.0, 0.0]
cylinder1Display.PolarAxes.Orientation = [90.0, 0.0, 0.0]
cylinder1Display.Position = [0.0, 0.0, -cylinderHalfHeight]
cylinder1Display.DataAxesGrid.Position = [0.0, 0.0, -cylinderHalfHeight]
cylinder1Display.PolarAxes.Translation = [0.0, 0.0, -cylinderHalfHeight]

renderView.Update()

theta = 18.0/180.0*np.pi
sin_theta = np.sin(theta)
numFrames = len(files)
t = np.linspace(0.0, 2*np.pi, numFrames, endpoint=False)
print(len(t), len(files))
gamma = -np.arcsin(-sin_theta*np.sin(t))
beta = -np.arcsin(sin_theta*np.cos(t)/np.sqrt(1.0 - (sin_theta*np.sin(t))**2))

# repeat sequence a few times
repetitions = 3
beta = np.tile(beta, repetitions)
gamma = np.tile(gamma, repetitions)

renderView.ResetCamera()
camera = renderView.GetActiveCamera()
camera.Elevation(-80)

Interact()

animationScene.PlayMode = 'Sequence'
animationScene.StartTime = 0
animationScene.NumberOfFrames = len(files)

animationScene.GoToFirst()
for b, g in zip(beta, gamma):
    transform2.Transform.Rotate = [g/np.pi*180, b/np.pi*180, 0.0]
    transformWater.Transform.Rotate = [g/np.pi*180, b/np.pi*180, 0.0]
    renderView.Update()
    time.sleep(0.1)
    curTime = animationScene.TimeKeeper.Time
    i = int(curTime)
    print(f"Time step: {i}")
    if i == len(files)-1:
        animationScene.GoToFirst()
    else:
        animationScene.GoToNext()
