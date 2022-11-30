#!/usr/bin/env pvpython

from paraview.simple import *
import numpy as np

# create a new 'Legacy VTK Reader'
path = '/Users/pumbaa/Sync/Data/LabOnAChip/8March22/'
filename = "chip.vtk" #"synth.vtk" # "chip_void.vtk"
chip_voidvtk = LegacyVTKReader(
    registrationName='chip_void',
    FileNames=[path + filename]
)

channel0 = LegacyVTKReader(
    registrationName='channel0',
    FileNames=[path + "channel_big.vtk"]
)

channel0s = LegacyVTKReader(
    registrationName='channel0s',
    FileNames=[path + "channel_small.vtk"]
)

data_base_path_big = "/Users/pumbaa/dune-master/dumux-microfluidic/build-cmake/appl/tilting-chip/big/"
data_base_path_small = "/Users/pumbaa/dune-master/dumux-microfluidic/build-cmake/appl/tilting-chip/small/"

transformCh00 = Transform(Input=channel0)
#transformCh00.Transform.Translate = [7.3, -2.0, 3.75]
transformCh00.Transform.Translate = [0.0, 0.0, 0.0]
transformCh00.Transform.Rotate = [0.0, 0.0, 0.0]

transformCh00s = Transform(Input=channel0s)
transformCh00s.Transform.Translate = [0.0, 0.0, 0.0]
transformCh00s.Transform.Rotate = [0.0, 0.0, 0.0]

transformCh0 = Transform(Input=transformCh00)
transformCh0.Transform.Rotate = [0.0, 0.0, 0.0]

transformCh0s = Transform(Input=transformCh00s)
transformCh0s.Transform.Rotate = [0.0, 0.0, 0.0]

transformCh11 = Transform(Input=channel0)
transformCh11.Transform.Translate = [0.0, 0.0, 0.0]
transformCh111 = Transform(Input=transformCh11)
transformCh111.Transform.Rotate = [0.0, 0.0, 180.0]

transformCh1 = Transform(Input=transformCh111)
transformCh1.Transform.Rotate = [0.0, 0.0, 0.0]

transformCh11s = Transform(Input=channel0s)
transformCh11s.Transform.Translate = [0.0, 0.0, 0.0]
transformCh111s = Transform(Input=transformCh11s)
transformCh111s.Transform.Rotate = [0.0, 0.0, 180.0]

transformCh1s = Transform(Input=transformCh111s)
transformCh1s.Transform.Rotate = [0.0, 0.0, 0.0]

steps = 239
files0 = [data_base_path_big + f"intersections-reservoir_0-{i}.vtu" for i in range(steps)]
files1 = [data_base_path_big + f"intersections-reservoir_1-{i}.vtu" for i in range(steps)]
intersections0 = XMLUnstructuredGridReader(registrationName='intersections0', FileName=files0)
intersections1 = XMLUnstructuredGridReader(registrationName='intersections1', FileName=files1)

files0s = [data_base_path_small + f"intersections-reservoir_0-{i}.vtu" for i in range(steps)]
files1s = [data_base_path_small + f"intersections-reservoir_1-{i}.vtu" for i in range(steps)]
intersections0s = XMLUnstructuredGridReader(registrationName='intersections0s', FileName=files0s)
intersections1s = XMLUnstructuredGridReader(registrationName='intersections1s', FileName=files1s)

animationScene = GetAnimationScene()
animationScene.UpdateAnimationUsingDataTimeSteps()

intersections0.TimeArray = 'None'
intersections1.TimeArray = 'None'
intersections0s.TimeArray = 'None'
intersections1s.TimeArray = 'None'

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
disp = GetDisplayProperties(transformCh0s, view=renderView)
disp.AmbientColor = [0.6666666666666666, 1.0, 1.0]
disp.DiffuseColor = [0.6666666666666666, 1.0, 1.0]
disp = GetDisplayProperties(transformCh1s, view=renderView)
disp.AmbientColor = [0.6666666666666666, 1.0, 1.0]
disp.DiffuseColor = [0.6666666666666666, 1.0, 1.0]

# translate to origin (somehow the geometry is not centered)
transform = Transform(Input=chip_voidvtk)
if "synth" in filename or filename == "chip.vtk":
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

transformWater0s = Transform(Input=intersections0s)
transformWater0s.Transform.Rotate = [0.0, 0.0, 0.0]

Show()
dispWater = GetDisplayProperties(transformWater0s, view=renderView)
dispWater.AmbientColor = [0.6666666666666666, 1.0, 1.0]
dispWater.DiffuseColor = [0.6666666666666666, 1.0, 1.0]

transformWater1 = Transform(Input=intersections1)
transformWater1.Transform.Rotate = [0.0, 0.0, 180.0]

Show()
dispWater = GetDisplayProperties(transformWater1, view=renderView)
dispWater.AmbientColor = [0.6666666666666666, 1.0, 1.0]
dispWater.DiffuseColor = [0.6666666666666666, 1.0, 1.0]

transformWater1s = Transform(Input=intersections1s)
transformWater1s.Transform.Rotate = [0.0, 0.0, 180.0]

Show()
dispWater = GetDisplayProperties(transformWater1s, view=renderView)
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
    with open(data_base_path_big + f"intersections-reservoir_0-{i}.txt") as metaData:
        g, b = metaData.readlines()[0].split()[:2]
        g, b = float(g), float(b)
        transformCh0.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transformCh1.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transformCh0s.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transformCh1s.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transform2.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transformWater0.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]
        transformWater0s.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 0.0]

    with open(data_base_path_big + f"intersections-reservoir_1-{i}.txt") as metaData:
        g, b = metaData.readlines()[0].split()[:2]
        g, b = float(g), float(b)
        transformWater1.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 180.0]
        transformWater1s.Transform.Rotate = [-g/np.pi*180, b/np.pi*180, 180.0]

    renderView.Update()

update(0)
Interact()

animationScene.PlayMode = 'Sequence'
animationScene.StartTime = 0
animationScene.NumberOfFrames = len(files0)

animationScene.GoToFirst()
for i in range(1, len(files0)):
    update(i)
    WriteImage(f"chip-anim_both-{i:04}.png")
    #time.sleep(0.05)
    curTime = animationScene.TimeKeeper.Time
    i = int(curTime)
    print(f"Time step: {i}")
    # if i == len(files0)-1:
    #     animationScene.GoToFirst()
    # else:
    animationScene.GoToNext()
