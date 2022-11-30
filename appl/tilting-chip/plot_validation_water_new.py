#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.gridspec as gridspec

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

data_filename = "../../../data/flow_new/piv-26-09-22.csv"
all_data = np.genfromtxt(data_filename, delimiter=",", dtype=np.float64, skip_header=1)
data = [
    all_data[:, [0, i]] for i in range(1, all_data.shape[1])
]

data_structured = {}
params_s_t = [(2, 18), (2, 15), (2, 10), (2, 7), (3, 18), (3, 15), (3, 10), (3, 7), (4, 18), (4, 15), (4, 10), (4, 7), (5, 18), (5, 15), (5, 10), (5, 7)]
varied_param = params_s_t

for d, labels in zip(data, varied_param):
    speed, angle = labels
    t = d[:,0]
    f = d[:,1]
    print(t, f)

    t = t[t < 60.0/float(speed)]
    f = f[:len(t)]

    # compute normalization using integral under curve which should be 300µl
    np.nan_to_num(f, copy=False)
    integral = np.trapz(f, dx=t[2]-t[1])
    data_structured[labels] = {}

    for vol in [150.0, 200.0, 250.0, 300.0]:
        factor = vol/integral
        f2 = f*factor
        t2 = t[t < 60.0/float(speed)]
        f2 = f2[:len(t2)]
        data_structured[labels][int(vol)] = (t2, f2)


def plot_measurement(ax, label, shift=4.0, everynth=1, vol=300):
    samples_per_second = len(data_structured[label][vol][0])/20.0
    shift = int(shift*samples_per_second)
    ax.plot(
        data_structured[label][vol][0][2:-2][::everynth],
        np.roll(running_mean(data_structured[label][vol][1], 5), shift)[::everynth],
        "-",  label=f"Measurement",
        ms=2*72.0/fig.dpi,
        alpha=0.5,
    )
    ax.set_xlim([0, 60.0/float(label[0])])


def plot_simulation_speed(ax, label, vol=300):
    #for crit_pressure in ["0", "10", "20", "30", "40"][::-1]:
    for crit_pressure in ["0","40"][::-1]:
    #for vol in [150, 200, 250, 300,]:
    #for vol in [250,]:
        speed, angle = label
        file_name = f"validation_s{speed}_a{angle}_p{crit_pressure}_d500_{vol}ul_t1_big.txt"

        label_str = f"Simulation-{crit_pressure}Pa"
        #label_str = f"Simulation"
        if Path(file_name).is_file():
            print(file_name)
            sim = np.genfromtxt(file_name, delimiter=" ", skip_header=1, skip_footer=2)
            print(sim[:,0][2]-sim[:,0][1])
            print(np.trapz(np.abs(sim[:,4]), dx=sim[:,0][2]-sim[:,0][1]))
            p0 = ax.plot(sim[:,0], -sim[:,4], "-", label=label_str, ms=72.0/fig.dpi)[0]
    ax.set_xlim([0, 60.0/float(speed)])
    ax.set_title(f"{speed} rpm, {angle}°, {vol}µl")
    #ax.plot(sim[:,0], sim[:,5], "-", label=label_str, ms=72.0/fig.dpi, color=p0.get_color(), alpha=0.2)


fig = plt.figure()
gs0 = gridspec.GridSpec(2, 1, figure=fig)
gs00 = gs0[0].subgridspec(1, 1)
gs01 = gs0[1].subgridspec(1, 2, width_ratios=[1, 7.0/9.0])

#fig, axes = plt.subplots(1, 2, figsize=(12,4), dpi=150, sharey=True, gridspec_kw={'width_ratios': [1, 4.0/7.0]})
#ax0, ax1, ax2, ax3 = axes
ax0 = fig.add_subplot(gs00[0])
ax1 = fig.add_subplot(gs01[0])
ax3 = fig.add_subplot(gs01[1], sharey=ax1)
plt.setp(ax3.get_yticklabels(), visible=False)
offset = 0.0
everynth = 10

plot_simulation_speed(ax3, (5, 18), vol=250)
plot_measurement(ax3, (5, 18), shift=6.5, everynth=everynth, vol=250)

plot_simulation_speed(ax0, (3, 18), vol=250)
plot_measurement(ax0, (3, 18), shift=15.7, everynth=everynth, vol=250)

plot_simulation_speed(ax1, (4, 18), vol=250)
plot_measurement(ax1, (4, 18), shift=13.0, everynth=everynth, vol=250)

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

set_size(4, 4, ax0)
set_size(4, 4, ax1)
#set_size(8, 4*0.5, ax2)
set_size(4, 4, ax3)

for ax in [ax0, ax1, ax3]:
    ax0.set_ylabel("Q in µl/s")
    ax1.set_ylabel("Q in µl/s")
    ax.set_xlabel("t in s")
    ax0.legend(frameon=False)
fig.tight_layout()
fig.savefig("speed_comparison_new.pdf", dpi=600)
plt.show()
