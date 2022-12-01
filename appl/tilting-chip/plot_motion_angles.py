import numpy as np
import matplotlib.pyplot as plt

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

t, tt, time, a, b = np.genfromtxt("angles.txt").T
tilt = np.arccos(np.cos(a)*np.cos(b))
fig, axes = plt.subplots(1, 2, figsize=(8,4))
axes[0].plot(t[:-1], tt[:-1], label="time warp function", color=colors[0])
axes[0].plot(t[:-1], t[:-1], "--", alpha=0.5, label="const. rotation velocity", color=colors[1])
axes[0].set_aspect(1.0)
axes[0].legend()
axes[1].plot(time, a*180/np.pi, label="pitch", color=colors[2])
axes[1].plot(time, b*180/np.pi, label="roll", color=colors[3])
axes[1].plot(time, tilt*180/np.pi, "--", alpha=0.5, label="tilt", color="k")
axes[1].set_xlabel("time in s")
axes[1].set_ylabel("angle in degree")
axes[1].legend()
fig.tight_layout()
plt.show()
