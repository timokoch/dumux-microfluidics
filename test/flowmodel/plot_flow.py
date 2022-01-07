import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(8,4))

def plot(filename, label):
    data = np.genfromtxt(filename, skip_header=1).T
    p0 = ax.plot(data[0,:], data[4,:], label=f"[ch0] {label}")[0]
    p1 = ax.plot(data[0,:], -data[5,:], "--", color=p0.get_color(), label=f"[ch1] {label}")
    #p0 = ax.plot(data[0,:], data[8,:], label=f"[ch0] {label}")[0]
    #p1 = ax.plot(data[0,:], -data[9,:], "--", color=p0.get_color(), label=f"[ch1] {label}")

plot("rpm2_10deg.txt", r"2rpm $\theta$=10°")
plot("rpm2_12deg.txt", r"2rpm $\theta$=12°")
plot("rpm2_14deg.txt", r"2rpm $\theta$=14°")
plot("rpm2_16deg.txt", r"2rpm $\theta$=16°")
plot("rpm2_18deg.txt", r"2rpm $\theta$=18°")
plot("rpm2_19deg.txt", r"2rpm $\theta$=19°")

# plot("rpm4_10deg.txt", r"4rpm $\theta$=10°")
# plot("rpm4_12deg.txt", r"4rpm $\theta$=12°")
# plot("rpm4_14deg.txt", r"4rpm $\theta$=14°")
# plot("rpm4_16deg.txt", r"4rpm $\theta$=16°")
# plot("rpm4_18deg.txt", r"4rpm $\theta$=18°")
# plot("rpm4_19deg.txt", r"4rpm $\theta$=19°")

#plot("rpm10_18deg.txt", r"10rpm $\theta$=18°")

ax.set_ylabel(r"Q in $\mu$l/s ${}^{**}$")
ax.set_xlabel(r"Time in s")
ax.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1.15), fancybox=False, frameon=False)

ax2 = ax.twinx()
mn, mx = ax.get_ylim()
ax2.set_ylim(mn/9e-1*0.026e-1, mx/9e-1*0.026e-1)
ax2.set_ylabel(r"$v_\mathrm{max}$ in m/s")

ax.set_title(
    r"$\mu$ = 1.0 $\cdot$10$^{-3}~$Pa$\cdot$s, "
    r"$R_{ch}$= 9.0$\cdot$10$^{-10}~\mu$l/(Pa s), "
    r"$V=300~\mu$l, $V_{A,0}$=266.53$~\mu$l, $V_{B,0}$=0.0$~\mu$l, $V_{ch0/1}$ = 16.735$~\mu$l",
    fontsize=8
)

plt.figtext(0.02, 0.02,
    r"${}^{**}$ negative Q means backflow",
    fontsize=8, horizontalalignment='left'
)

plt.tight_layout()
fig.savefig("flow_chip_rpm2.pdf", dpi=300)

plt.show()
