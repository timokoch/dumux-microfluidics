#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def roll(vec, gamma):
    s, c, one = np.sin(gamma), np.cos(gamma), np.ones(gamma.shape)
    x = one*vec[0]
    y = c*vec[1] - s*vec[2]
    z = s*vec[1] + c*vec[2]
    return np.array([x, y, z])

def pitch(vec, beta):
    s, c, one = np.sin(beta), np.cos(beta), np.ones(beta.shape)
    x = c*vec[0] + s*vec[2]
    y = one*vec[1]
    z = -s*vec[0] + c*vec[2]
    return np.array([x, y, z])

def tilt_angle(vec):
    """
        The angle between normal vector of platform the z-axis unit vector
        That is via scalar product: cos(angle) = vec*ez
    """
    return np.arccos(vec[2])

def compute_tilt(beta, gamma):
    """
        For current tilt: roll around x-axis by gamma, pitch around
        y-axis by beta and then compute angle to z-axis
    """
    return tilt_angle(pitch(roll(np.array([0, 0, 1]), gamma), beta))

def compute_tilt_simplified(beta, gamma):
    """
        Simplified but exact formula
    """
    return np.arccos(np.cos(beta)*np.cos(gamma))

angle = 19.0
theta = angle/180.0*np.pi
sinTheta = np.sin(theta)
rpm = 5.23

time = np.linspace(0, 60.0/rpm, 500)
t = time*2.0*np.pi*rpm/60.0

sint = np.sin(t)
cost = np.cos(t)

# this gives perfectly constant tilt over time
gamma = np.arcsin(-sinTheta*sint)
beta = np.arcsin(-sinTheta*cost/(np.sqrt(1.0 - sinTheta**2*sint**2)))

# this is an approximation for small angles
gammaApprox = -theta*sint
betaApprox = -theta*cost

tilt = compute_tilt(beta, gamma)
tilt2 = compute_tilt_simplified(beta, gamma)
tiltApprox = compute_tilt(betaApprox, gammaApprox)

# plot
fig, (ax0, axt, ax1) = plt.subplots(3, 1, figsize=(10,10), sharex=True)

# approximations
ax0.plot(time[::5], betaApprox[::5]*180.0/np.pi, "+", label="approx beta", )
ax0.plot(time[::5], gammaApprox[::5]*180.0/np.pi, "+", label="approx gamma")

axt.plot(time[::5], tiltApprox[::5]*180.0/np.pi, "+", label="tilt approx")

# exact (this gives perfectly constant tilt over time)
ax0.plot(time, beta*180.0/np.pi, label="beta")
ax0.plot(time, gamma*180.0/np.pi, label="gamma")

axt.plot(time, tilt*180.0/np.pi, label="tilt")
axt.plot(time, tilt2*180.0/np.pi, label="tilt2")

# measurements
mt, mx, my, mz = measured_data = np.genfromtxt("gyroscope.csv", delimiter=",", skip_header=1, unpack=True)
mt = mt/1000.0 # to seconds
# mt = mt - mt[0]
# freqs = np.fft.fftfreq(len(w), d=mt)
# print(freqs.min(), freqs.max())

# # Find the peak in the coefficients
# idx = np.argmax(np.abs(w))
# freq = freqs[idx]
# freq_in_hertz = np.abs(freq * len(mt)/mt[-1])
# print(freq_in_hertz)

def mask_and_shift(mt, mx, my, mz, start_time=5.0):
    period = 60.0/rpm
    mask_m = (mt >= start_time) & (mt < start_time+period) # pick one period
    mt, mx, my, mz = mt[mask_m], mx[mask_m], my[mask_m], mz[mask_m]
    shift = -np.argmin(my)
    mx, my, mz = np.roll(mx, shift), np.roll(my, shift), np.roll(mz, shift)
    mt = mt - mt[0]
    return mt, mx, my, mz

mt, mx, my, mz = mask_and_shift(mt, mx, my, mz, 8.1)

ax0.plot(mt[::3], mz[::3], "kd", label="gyroscope gamma", markersize=2)
ax0.plot(mt[::3], my[::3], "ks", label="gyroscope beta", markersize=2)

axt.plot(mt[::3], compute_tilt_simplified(my[::3]*np.pi/180.0, mz[::3]*np.pi/180.0)*180.0/np.pi, "k+", label="gyroscope tilt", markersize=2)

# errors
ax1.plot(time, (beta-betaApprox)/np.max(beta)*100, label="beta")
ax1.plot(time, (gamma-gammaApprox)/np.max(gamma)*100, label="gamma")
ax1.plot(time, (tilt-tiltApprox)/np.max(tilt)*100, label="tilt")

ax0.set_ylabel("Angle in degree")
axt.set_ylabel("Angle in degree")
ax1.set_ylabel("Relative error of approximation in %")
ax1.set_xlabel("Time in seconds")
ax0.legend()
axt.legend()
ax1.legend()

ax0.set_title(f"{rpm} rpm and {angle}° tilt")

fig.tight_layout()
fig.savefig("tilting_angles.pdf")

plt.show()

