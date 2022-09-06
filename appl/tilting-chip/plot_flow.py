import matplotlib.pyplot as plt
import numpy as np
import sys

# plotting script for debugging/visualizing a result
filename = str(sys.argv[1])

fig, (ax, ax3) = plt.subplots(2, 1, figsize=(8,8))

def plot(filename):
    data = np.genfromtxt(filename, skip_header=1).T
    p0 = ax.plot(data[0, :], -data[4, :], "-", label=f"[Channel 0]")[0]
    ax.plot(data[0,:], data[5,:], "--", color=p0.get_color(), label=f"[Channel 1]", alpha=0.2)
    ax3.plot(data[0,:], data[2,:], "-", label=f"[Reservoir 0]")
    ax3.plot(data[0,:], data[3,:], "-", label=f"[Reservoir 1]")

plot(filename)

ax.set_title(filename)
ax.set_ylabel(r"Q in $\mu$l/s")
ax.set_xlabel(r"Time in s")
ax3.set_ylabel(r"Volume in $\mu$l")
ax3.set_xlabel(r"Time in s")
ax.legend(ncol=2, loc='lower center', bbox_to_anchor=(0.5, 1.15), fancybox=False, frameon=False)
ax3.legend(ncol=2, loc='lower center', bbox_to_anchor=(0.5, 1.15), fancybox=False, frameon=False)

# plt.figtext(0.02, 0.02,
#     r"${}^{**}$ negative Q or $v_\mathrm{max}$ means backflow",
#     fontsize=8, horizontalalignment='left'
# )

plt.tight_layout()
fig.savefig("test.pdf", dpi=300)

plt.show()
