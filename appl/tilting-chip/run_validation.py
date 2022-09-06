#!/usr/bin/env python3

import subprocess
import params
import itertools
import multiprocessing as mp
from params import make_cmd

# we run for different speeds and critical pressures
# these are the tilting speeds we have data for
tilting_speeds = ["3", "5", "8"]
# we test for different critical pressures
# to investigate the influence of the stop valve effect
#crit_pressures = ["0", "10", "20", "40"]
crit_pressures = ["30",]
t_factor = [2,]
vol = [200, 250, 300,]

cases = [
    tilting_speeds,
    crit_pressures,
    t_factor,
    vol,
]

def run_case(case):
    speed, crit_pressure, tf, vol = case

    args = params.DEFAULT_CASE_VALIDATION
    args["Problem.CapillaryStopPressure"] = crit_pressure
    dryingT = 500.0 # micron
    args["Problem.DryingThresholdInMilliMeter"] = f"{dryingT/1000.0}"
    args["Problem.RotationsPerMinute"] = speed
    args["Problem.ChannelTransmissibility"] = str(float(tf)*float(args["Problem.ChannelTransmissibility"]))
    args["Problem.InitialVolumeInMicroLiter"] = str(vol)

    cmd = make_cmd(args, f"validation_s{speed}_p{crit_pressure}_d{int(dryingT)}_{vol}ul_t{tf}_big.txt")
    subprocess.run(cmd)

if __name__ == '__main__':
    with mp.Pool(8) as p:
        p.map(run_case, itertools.product(*cases))
