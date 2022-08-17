import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(8,4))

def plot(filename, lab):
    data = np.genfromtxt(filename, skip_header=1).T
    p0 = ax.plot(data[0, :], -data[4, :], "*", label=f"[ch0] {lab}")[0]
    ax.plot(data[0,:], data[5,:], "--", color=p0.get_color(), label=f"[ch1] {lab}", alpha=0.2)

plot("18_5_500_big", r"5rpm $\theta$=18° 500µl (outer)")
plot("18_5_400_big", r"5rpm $\theta$=18° 400µl (outer)")
plot("18_5_300_big", r"5rpm $\theta$=18° 300µl (outer)")
#plot("18_5_200_big", r"5rpm $\theta$=18° 200µl (outer)")

#plot("18_4_500_big", r"4rpm $\theta$=18° 500µl (outer)")
#plot("18_4_400_big", r"4rpm $\theta$=18° 400µl (outer)")
#plot("18_4_300_big", r"4rpm $\theta$=18° 300µl (outer)")
#plot("18_4_200_big", r"4rpm $\theta$=18° 200µl (inner)")

#plot("18_5_200_big", r"5rpm $\theta$=18° 200µl (outer)")
#plot("18_5_250_small", r"5rpm $\theta$=18° 250µl (inner)")
#plot("18_5_200_big", r"5rpm $\theta$=18° 200µl (outer)")
#plot("18_5_200_small", r"5rpm $\theta$=18° 200µl (inner)")

ax.set_ylabel(r"Q in $\mu$l/s ${}^{**}$")
ax.set_xlabel(r"Time in s")
ax.legend(ncol=2, loc='lower center', bbox_to_anchor=(0.5, 1.15), fancybox=False, frameon=False)
ax2 = ax.twinx()

mn, mx = ax.get_ylim()
# there is a constant factor converting between max velocity and flow rate
ax2.set_ylim(mn*1e-9*1117183.943187024, mx*1e-9*1117183.943187024)
ax2.set_ylabel(r"$v_\mathrm{max}$ in m/s")

plt.figtext(0.02, 0.02,
    r"${}^{**}$ negative Q or $v_\mathrm{max}$ means backflow",
    fontsize=8, horizontalalignment='left'
)

plt.tight_layout()
fig.savefig("new_chip_rpm5_both.pdf", dpi=300)

plt.show()
