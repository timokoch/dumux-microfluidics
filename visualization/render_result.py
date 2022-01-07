#!/usr/bin/env pvpython

from paraview.simple import *
import numpy as np

# create a new 'Legacy VTK Reader'
filename = "chip_void.vtk" #"synth.vtk" # "chip_void.vtk"
chip_voidvtk = LegacyVTKReader(
    registrationName='chip_void',
    FileNames=['/Users/pumbaa/Sync/Data/LabOnAChip/' + filename]
)

channel0 = LegacyVTKReader(
    registrationName='channel0',
    FileNames=["channel.vtk"]
)

transformCh00 = Transform(Input=channel0)
transformCh00.Transform.Translate = [7.3, -2.0, 3.75]
transformCh00.Transform.Rotate = [0.0, 0.0, 0.0]

transformCh0 = Transform(Input=transformCh00)
transformCh0.Transform.Rotate = [0.0, 0.0, 0.0]

transformCh11 = Transform(Input=channel0)
transformCh11.Transform.Translate = [7.3, -2.0, 3.75]
transformCh111 = Transform(Input=transformCh11)
transformCh111.Transform.Rotate = [0.0, 0.0, 180.0]

transformCh1 = Transform(Input=transformCh111)
transformCh1.Transform.Rotate = [0.0, 0.0, 0.0]

steps = 200
files0 = [f"/Users/pumbaa/dune-master/dumux-microfluidic/build-cmake/test/flowmodel/intersections-reservoir_0-{i}.vtu" for i in range(steps)]
files1 = [f"/Users/pumbaa/dune-master/dumux-microfluidic/build-cmake/test/flowmodel/intersections-reservoir_1-{i}.vtu" for i in range(steps)]
intersections0 = XMLUnstructuredGridReader(registrationName='intersections0', FileName=files0)
intersections1 = XMLUnstructuredGridReader(registrationName='intersections1', FileName=files1)

animationScene = GetAnimationScene()
animationScene.UpdateAnimationUsingDataTimeSteps()

intersections0.TimeArray = 'None'
intersections1.TimeArray = 'None'

# get active view
renderView = GetActiveViewOrCreate('RenderView')
renderView.Background = [1.0, 1.0, 1.0]

# make channel blue
disp = GetDisplayProperties(transformCh0, view=renderView)
disp.AmbientColor = [0.6666666666666666, 1.0, 1.0]
disp.DiffuseColor = [0.6666666666666666, 1.0, 1.0]
disp = GetDisplayProperties(transformCh1, view=renderView)
disp.AmbientColor = [0.6666666666666666, 1.0, 1.0]
disp.DiffuseColor = [0.6666666666666666, 1.0, 1.0]

# translate to origin (somehow the geometry is not centered)
transform = Transform(Input=chip_voidvtk)
if "synth" in filename:
    transform.Transform.Translate = [0.0, 0.0, 0.0]
else:
    transform.Transform.Translate = [7.3, -2.0, 3.75]

transform2 = Transform(Input=transform)
transform2.Transform.Rotate = [0.0, 0.0, 0.0]

Show()
disp = GetDisplayProperties(transform2, view=renderView)
disp.Opacity = 0.15

transformWater0 = Transform(Input=intersections0)
transformWater0.Transform.Rotate = [0.0, 0.0, 0.0]

Show()
dispWater = GetDisplayProperties(transformWater0, view=renderView)
dispWater.AmbientColor = [0.6666666666666666, 1.0, 1.0]
dispWater.DiffuseColor = [0.6666666666666666, 1.0, 1.0]

transformWater1 = Transform(Input=intersections1)
transformWater1.Transform.Rotate = [0.0, 0.0, 180.0]

Show()
dispWater = GetDisplayProperties(transformWater1, view=renderView)
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

renderView.ResetCamera()
camera = renderView.GetActiveCamera()
camera.Elevation(-80)

def update(i):
    with open(f"/Users/pumbaa/dune-master/dumux-microfluidic/build-cmake/test/flowmodel/intersections-reservoir_0-{i}.txt") as metaData:
        g, b = metaData.readlines()[0].split()[:2]
        g, b = float(g), float(b)
        transformCh0.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transformCh1.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transform2.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transformWater0.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]

    with open(f"/Users/pumbaa/dune-master/dumux-microfluidic/build-cmake/test/flowmodel/intersections-reservoir_1-{i}.txt") as metaData:
        g, b = metaData.readlines()[0].split()[:2]
        g, b = float(g), float(b)
        transformWater1.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 180.0]

    renderView.Update()

update(0)
Interact()

animationScene.PlayMode = 'Sequence'
animationScene.StartTime = 0
animationScene.NumberOfFrames = len(files0)

animationScene.GoToFirst()
for i in range(1, len(files0)):
    update(i)
    WriteImage(f"chip-anim-{i:04}.png")
    #time.sleep(0.05)
    curTime = animationScene.TimeKeeper.Time
    i = int(curTime)
    print(f"Time step: {i}")
    # if i == len(files0)-1:
    #     animationScene.GoToFirst()
    # else:
    animationScene.GoToNext()
