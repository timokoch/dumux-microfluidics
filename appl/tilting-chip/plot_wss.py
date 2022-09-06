#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

data = [
    np.genfromtxt("../../../data/flow/3rpm_300ul_t19_cellmedia.csv", delimiter=","),
    np.genfromtxt("../../../data/flow/5rpm_300ul_t19_cellmedia.csv", delimiter=","),
    np.genfromtxt("../../../data/flow/8rpm_300ul_t19_cellmedia.csv", delimiter=","),
]
data_structured = {}
speed = [3, 5, 8]
for d, label in zip(data, speed):
    t = d[:,0]
    f = d[:,1]
    # compute normalization using integral under curve which should be 300µl
    np.nan_to_num(f, copy=False)
    integral = np.trapz(f, dx=20.0/len(f))
    factor = (300.0)/integral
    print(factor)
    f = f*factor
    t = t[t < 60.0/float(label)]
    f = f[:len(t)]
    data_structured[f"{label}"] = (t, f)


def plot_measurement(ax, label, shift=4.0, everynth=1):
    samples_per_second = len(data_structured[label][0])/20.0
    shift = int(shift*samples_per_second)
    ax.plot(
        data_structured[label][0][2:-2][::everynth],
        np.roll(running_mean(data_structured[label][1], 5), shift)[::everynth],
        ".",  label=f"M-{label}rpm",
        ms=2*72.0/fig.dpi,
        alpha=0.5,
    )
    ax.set_xlim([0, 60.0/float(label)])


volume = ["100", "150", "200", "250", "300", "400", "800"]
def plot_simulation_speed(ax, label, mod=False):
    for v in volume:
        file_name = f"wss_s{label}_p0_{v}ul_d500_300_big.txt"
        print(f"Parse {file_name}")
        label_str= f"S-{label}rpm|CV" if mod else f"S-{label}rpm-{v}uL"
        if Path(file_name).is_file():
            sim = np.genfromtxt(file_name, delimiter=" ", skip_header=1, skip_footer=2)
            ax.plot(sim[:,0], -sim[:,4], "-", label=label_str, ms=72.0/fig.dpi)[0]
        else:
            print(f"-- File {file_name} not found")
    ax.set_xlim([0, 60.0/float(label)])


fig, axes = plt.subplots(1, 3, figsize=(12,4), dpi=150, sharey=True)
ax0, ax1, ax2 = axes

plot_simulation_speed(ax0, "3")
plot_simulation_speed(ax1, "5")
plot_simulation_speed(ax2, "8")

offset = 2.0
everynth = 10
plot_measurement(ax0, "3", shift=18.0+offset, everynth=everynth)
plot_measurement(ax1, "5", shift=19.5+offset, everynth=everynth)
plot_measurement(ax2, "8", shift=11.0+offset, everynth=everynth)

for ax in axes.flatten():
    ax.set_ylabel(r"WSS in dyne/cm${}^2$")
    ax.set_xlabel("t in s")
    ax.legend()
fig.tight_layout()
fig.savefig("wss_volume_comparison.pdf", dpi=600)
plt.show()
