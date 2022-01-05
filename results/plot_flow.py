import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(8,4))

def plot(filename, label):
    data = np.genfromtxt(filename, skip_header=1).T
    p0 = ax.plot(data[0,:], data[3,:], label=f"[{label}] channel 0")[0]
    p1 = ax.plot(data[0,:], -data[4,:], "--", color=p0.get_color(), label=f"[{label}] channel 1")

plot("output_03.txt", r"0.3 s$^{-1}$")
plot("output_03B.txt", r"0.3 s$^{-1}$ R$\downarrow$")
plot("output_03C.txt", r"0.3 s$^{-1}$ V$\downarrow$")
plot("output_04.txt", r"0.4 s$^{-1}$")
plot("output_05.txt", r"0.5 s$^{-1}$")

ax.set_ylabel(r"Q in mm$^3$/s ${}^{**}$")
ax.set_xlabel(r"Time in s")

plt.figtext(0.01, 0.01, r"${}^{**}$ negative Q means backflow", fontsize=8, horizontalalignment='left')
plt.legend(ncol=3)
plt.tight_layout()
fig.savefig("flow_chip.pdf", dpi=300)

plt.show()
