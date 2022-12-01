import numpy as np
import matplotlib.pyplot as plt

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

t, tt, time, a, b = np.genfromtxt("angles.txt").T
tilt = np.arccos(np.cos(a)*np.cos(b))
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
axes[0].plot(t[:-1], tt[:-1], label="time warp function", color=colors[0])
axes[0].plot(t[:-1], t[:-1], "--", alpha=0.5, label="const. rotation velocity", color=colors[1])
axes[0].set_aspect(1.0)
axes[0].legend()
axes[1].plot(t[:-1], (tt[1:]-tt[:-1])/(t[1:]-t[:-1]), label="relative velocity", color=colors[0])
axes[1].plot(t[:-1], np.ones(t[:-1].shape), "--", alpha=0.5, label="const. rotation velocity", color=colors[1])
axes[1].set_ylabel("relative speed-up")
axes[1].set_ylim([axes[1].get_ylim()[0], axes[1].get_ylim()[1]*1.2])
axes[1].legend()
axes[2].plot(time, a*180/np.pi, label="pitch", color=colors[2])
axes[2].plot(time, b*180/np.pi, label="roll", color=colors[3])
axes[2].plot(time, tilt*180/np.pi, "--", alpha=0.5, label="tilt", color="k")
axes[2].set_xlabel("time in s")
axes[2].set_ylabel("angle in degree")
axes[2].legend()
fig.tight_layout()
plt.show()
