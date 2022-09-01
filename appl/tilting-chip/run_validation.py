#!/usr/bin/env python3

import subprocess
import params
import multiprocessing as mp
from params import make_cmd

# we run for different speeds and critical pressures
# these are the tilting speeds we have data for
tilting_speeds = ["3", "5", "8"]
# we test for different critical pressures
# to investigate the influence of the stop valve effect
crit_pressures = ["20", "25", "30", "35", "40"]

def run_validation_case(crit_pressure):
    args = params.DEFAULT_CASE_VALIDATION
    args["Problem.CapillaryStopPressure"] = crit_pressure

    for speed in tilting_speeds:
        args["Problem.RotationsPerMinute"] = speed
        cmd = make_cmd(args, f"s{speed}_p{crit_pressure}_300_big.txt")
        subprocess.run(cmd)

if __name__ == '__main__':
    with mp.Pool(len(crit_pressures)) as p:
        p.map(run_validation_case, crit_pressures)
