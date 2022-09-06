#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

data = [
    np.genfromtxt("../../../data/flow/3rpm_300ul_t19_cellmedia.csv", delimiter=",", dtype=np.float64),
    np.genfromtxt("../../../data/flow/5rpm_300ul_t19_cellmedia.csv", delimiter=",", dtype=np.float64),
    np.genfromtxt("../../../data/flow/8rpm_300ul_t19_cellmedia.csv", delimiter=",", dtype=np.float64),
]
data_structured = {}
speed = [3, 5, 8]
for d, label in zip(data, speed):
    t = d[:,0]
    f = d[:,1]

    # compute normalization using integral under curve which should be 300µl
    np.nan_to_num(f, copy=False)
    integral = np.trapz(f, dx=t[2]-t[1])
    data_structured[f"{label}"] = {}

    for vol in [150.0, 200.0, 250.0, 300.0]:
        factor = vol/integral
        f2 = f*factor
        t2 = t[t < 60.0/float(label)]
        f2 = f2[:len(t2)]
        data_structured[f"{label}"][int(vol)] = (t2, f2)


def plot_measurement(ax, label, shift=4.0, everynth=1, vol=300):
    samples_per_second = len(data_structured[label][vol][0])/20.0
    shift = int(shift*samples_per_second)
    ax.plot(
        data_structured[label][vol][0][2:-2][::everynth],
        np.roll(running_mean(data_structured[label][vol][1], 5), shift)[::everynth],
        "-",  label=f"Meas",
        ms=2*72.0/fig.dpi,
        alpha=0.5,
    )
    ax.set_xlim([0, 60.0/float(label)])


def plot_simulation_speed(ax, label):
    #for crit_pressure in ["0", "10", "20", "30", "40"][::-1]:
    for crit_pressure in ["40",][::-1]:
    #for vol in [150, 200, 250, 300,]:
    #for vol in [250,]:
        vol = 250
        file_name = f"validation_s{label}_p{crit_pressure}_d700_{vol}ul_t1_big.txt"
        #label_str = f"S-{crit_pressure}Pa"
        label_str = f"Sim"
        if Path(file_name).is_file():
            print(file_name)
            sim = np.genfromtxt(file_name, delimiter=" ", skip_header=1, skip_footer=2)
            print(sim[:,0][2]-sim[:,0][1])
            print(np.trapz(np.abs(sim[:,4]), dx=sim[:,0][2]-sim[:,0][1]))
            p0 = ax.plot(sim[:,0], -sim[:,4], "-", label=label_str, ms=72.0/fig.dpi)[0]
    ax.set_xlim([0, 60.0/float(label)])
    ax.set_title(f"{label} rpm, {vol}µl")
    #ax.plot(sim[:,0], sim[:,5], "-", label=label_str, ms=72.0/fig.dpi, color=p0.get_color(), alpha=0.2)


fig, axes = plt.subplots(1, 2, figsize=(12,4), dpi=150, sharey=True, gridspec_kw={'width_ratios': [1, 3.0/8.0]})
#ax0, ax1, ax2, ax3 = axes
ax0, ax2 = axes
offset = 2.0
everynth = 10

plot_simulation_speed(ax0, "3")
#plot_measurement(ax0, "3", shift=18.0+offset, everynth=everynth, vol=300)
plot_measurement(ax0, "3", shift=18.0+offset, everynth=everynth, vol=250)
#plot_measurement(ax0, "3", shift=18.0+offset, everynth=everynth, vol=200)
#plot_measurement(ax0, "3", shift=18.0+offset, everynth=everynth, vol=150)

# plot_simulation_speed(ax1, "5")
# #plot_measurement(ax1, "5", shift=19.5+offset, everynth=everynth, vol=300)
# plot_measurement(ax1, "5", shift=19.5+offset, everynth=everynth, vol=250)
# #plot_measurement(ax1, "5", shift=19.5+offset, everynth=everynth, vol=200)
# #plot_measurement(ax1, "5", shift=19.5+offset, everynth=everynth, vol=150)

plot_simulation_speed(ax2, "8")
#plot_measurement(ax2, "8", shift=11.0+offset, everynth=everynth, vol=300)
plot_measurement(ax2, "8", shift=11.7+offset, everynth=everynth, vol=250)
#plot_measurement(ax2, "8", shift=11.0+offset, everynth=everynth, vol=200)
#plot_measurement(ax2, "8", shift=11.0+offset, everynth=everynth, vol=150)

# plot_simulation_speed(ax3, "3")
# plot_measurement(ax3, "3", shift=18.0+offset, everynth=everynth)

# plot_simulation_speed(ax3, "5")
# plot_measurement(ax3, "5", shift=19.5+offset, everynth=everynth)

# plot_simulation_speed(ax3, "8")
# plot_measurement(ax3, "8", shift=11.0+offset, everynth=everynth)

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

set_size(6, 6, ax0)
set_size(6, 6*3.0/8.0, ax2)

for ax in axes.flatten():
    ax0.set_ylabel("Q in µl/s")
    ax.set_xlabel("t in s")
    ax0.legend(frameon=False)
fig.tight_layout()
fig.savefig("speed_comparison.pdf", dpi=600)
plt.show()
