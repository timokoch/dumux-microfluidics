#!/usr/bin/env python3

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from pathlib import Path

tilting_speeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
reservoir = ["small", "big"]
#volume = np.linspace(100, 800, endpoint=True, num=10)
angles = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

X, Y = np.meshgrid(tilting_speeds, angles)

def max_wss(speed, tilt):
    def generate():
        for s, t in zip(speed, tilt):
            file_name = f"wss_s{s}_t{t}_p0_300ul_d500_300_big.txt"
            #file_name = f"wss_s{s}_t{t}_p0_200ul_d500_300_small.txt"
            print(f"Parse {file_name}")
            if Path(file_name).is_file():
                sim = np.genfromtxt(file_name, delimiter=" ", skip_header=1, skip_footer=2)
                yield np.max(sim[:,6])*10 # Pa to dyne/cm^2
                #yield np.max(-sim[:,4])
            else:
                yield 0
    return np.array(list(generate()))

Z = max_wss(X.flatten(), Y.flatten()).reshape(X.shape)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('Maximum WSS (outer circuit | 300µl)')
ax.set_xlabel('tilting speed (rpm)')
ax.set_ylabel('tilt angle (degree)')
ax.set_zlabel(r'Maximum WSS (dyne/cm$^{2}$)')
fig.tight_layout()
plt.show()
fig.savefig("wss_wireframe.pdf", dpi=600)
plt.show()
