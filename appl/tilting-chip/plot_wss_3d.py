#!/usr/bin/env python3

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib

tilting_speeds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#volume = np.linspace(100, 800, endpoint=True, num=10)
angles = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

#reservoir = "small"
reservoir = "big"

X, Y = np.meshgrid(tilting_speeds, angles)

def max_wss(speed, tilt):
    def generate():
        for s, t in zip(speed, tilt):
            if reservoir == "big":
                file_name = f"wss_s{s}_t{t}_p0_300ul_d500_big.txt"
            else:
                file_name = f"wss_s{s}_t{t}_p0_200ul_d500_small.txt"
            print(f"Parse {file_name}")
            if Path(file_name).is_file():
                sim = np.genfromtxt(file_name, delimiter=" ", skip_header=1, skip_footer=2)
                yield np.max(sim[:,6])*10 # Pa to dyne/cm^2
            else:
                yield 0
    return np.array(list(generate()))

def directionality(speed, tilt):
    def generate():
        for s, t in zip(speed, tilt):
            if reservoir == "big":
                file_name = f"wss_s{s}_t{t}_p0_300ul_d500_big.txt"
            else:
                file_name = f"wss_s{s}_t{t}_p0_200ul_d500_small.txt"
            print(f"Parse {file_name}")
            if Path(file_name).is_file():
                sim = np.genfromtxt(file_name, delimiter=" ", skip_header=1, skip_footer=2)
                yield np.trapz(np.abs(sim[:,4]), dx=sim[:,0][2]-sim[:,0][1])/300.0 # Pa to dyne/cm^2
            else:
                yield 0
    return np.array(list(generate()))

Z = max_wss(X.flatten(), Y.flatten()).reshape(X.shape)
Z2 = directionality(X.flatten(), Y.flatten()).reshape(X.shape)

color_dimension = Z2 # change to desired fourth dimension
minn, maxx = color_dimension.min(), color_dimension.max()
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='viridis')
m.set_array([])
fcolors = m.to_rgba(color_dimension)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                edgecolor='none', facecolors=fcolors, vmin=minn, vmax=maxx)

if reservoir == "big":
    ax.set_title('Maximum WSS (outer circuit / V=300µl)')
else:
    ax.set_title('Maximum WSS (inner circuit / V=200µl)')
ax.set_xlabel('tilting speed (rpm)')
ax.set_ylabel('tilt angle (degree)')
ax.set_zlabel(r'Maximum WSS (dyne/cm$^{2}$)')

if reservoir == "big":
    ax.set_zlim(0, 12)
else:
    ax.set_zlim(0, 8)
ax.view_init(elev=20, azim=-142)
cbar = fig.colorbar(surf, shrink=0.5, aspect=10)
cbar.set_label(r"$V\,/\,\int \! Q \,dt$")
fig.tight_layout()
plt.show()
fig.savefig("wss_wireframe.pdf", dpi=600)
plt.show()
