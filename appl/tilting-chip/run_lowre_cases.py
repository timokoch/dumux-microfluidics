#!/usr/bin/env python3

import subprocess
import params
import itertools
import multiprocessing as mp
import numpy as np
from params import make_cmd

# we run for different speeds and filling volumes and both reservoir types
tilting_speeds = [0.2,]
#reservoir = ["small", "big"]
#volume = np.linspace(100, 800, endpoint=True, num=10)
angles = [20,]

#tilting_speeds = [1,]
reservoir = ["big",]
#volume = np.linspace(100, 800, endpoint=True, num=10)
#angles = [10,]

cases = [
    tilting_speeds,
    reservoir,
    angles,
]

def run_case(case):
    speed, reservoir, angle = case
    args = params.DEFAULT_CASE_LOWRE_BIG if reservoir == "big" else params.DEFAULT_CASE_LOWRE_SMALL

    if reservoir == "big":
        args["Problem.InitialVolumeInMicroLiter"] = "300"
    else:
        args["Problem.InitialVolumeInMicroLiter"] = "200"

    args["Problem.MotionPointsX"] = "0.0 0.03 0.47 0.5"
    args["Problem.MotionPointsY"] = "0.0 0.18 0.43 0.5"

    args["Problem.CapillaryStopPressure"] = "0"
    drying_thres = 500.0 # micron
    args["Problem.DryingThresholdInMilliMeter"] = f"{drying_thres/1000.0}"
    args["Problem.RotationsPerMinute"] = speed
    args["TimeLoop.Dt"] = 1.0/(speed/60.0)/400.0
    args["Problem.Angle"] = angle

    cmd = make_cmd(args, f"lowre_s{speed}_t{angle}_p0_{args['Problem.InitialVolumeInMicroLiter']}ul_d{int(drying_thres)}_{reservoir}.txt")
    subprocess.run(cmd)

if __name__ == '__main__':
    with mp.Pool(64) as p:
        p.map(run_case, itertools.product(*cases))
