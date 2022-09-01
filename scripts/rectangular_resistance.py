#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

"""
Laminar flow rectangular tube (H*L):
see PDF in docs
"""

def plot(aspect_ratio, channel_length, scale, img, mu=1e-3):
    L = 1.0
    H = L/aspect_ratio
    G = 5000 #-dp/dz

    def dimless_velocity(y, z, num):
        h = 0.5*H
        w = 0.5*L
        x = y-w
        y = z-h
        yhat = y/h
        result = 1.0-yhat**2
        for n in range(1, num+1):
            beta = (2.0*n - 1.0)*np.pi/2
            prefactor = 4.0*((-1)**n)/beta**3/np.cosh(beta*aspect_ratio)
            result += prefactor*np.cosh(beta*x/h)*np.cos(beta*yhat)
        return result


    def velocity(y, z, num):
        HStar = H*scale
        return dimless_velocity(y, z, num)*(G*HStar**2)/(8*mu)


    def dimless_wall_shear_stress(y, z, num):
        h = 0.5*H
        w = 0.5*L
        x = y-w
        y = z-h
        yhat = y/h
        result = 0.0
        for n in range(1, num+1):
            beta = (2.0*n - 1.0)*np.pi/2
            prefactor = ((-1)**n)/beta**2/np.cosh(beta*aspect_ratio)
            result += prefactor*np.sinh(beta*x/h)*np.cos(beta*yhat)
        return result

    def wall_shear_stress(y, z, num):
        HStar = H*scale
        return dimless_wall_shear_stress(y, z, num)*(G*HStar)


    def dimless_flow_rate(num):
        result = 1.0
        factor = -6.0/aspect_ratio
        for n in range(1, num+1):
            beta = (2.0*n - 1.0)*np.pi/2.0
            prefactor = 1.0/(beta**5)
            result += factor*prefactor*np.tanh(beta*aspect_ratio)
        return result


    def flow_rate(num):
        LStar = L*scale
        HStar = H*scale
        return dimless_flow_rate(num)*(G*LStar*HStar**3)/(12*mu)


    def transmissibility(num, length):
        LStar = L*scale
        HStar = H*scale
        return dimless_flow_rate(num)*(LStar*HStar**3)/(12*mu*length)


    print("Channel volume (µl): ", H*L*scale**2*channel_length*1e9)
    print("Channel H (mm): ", H*scale*1e3)
    print("Channel W (mm): ", L*scale*1e3)
    print("Channel L (mm): ", channel_length*1e3)
    print("Pressure gradient: ", G)
    Q = flow_rate(10)
    print("Flow rate: ", Q)
    print("Transmissibility: ", transmissibility(10, length=channel_length))
    print("Transmissibility: ", transmissibility(10, length=channel_length)*channel_length)

    print("Max. velocity: ", velocity(0.5*L, 0.5*H, num=10))
    print("Max. Wall shear stress: ", wall_shear_stress(0, 0.5*H, num=10))
    print(f"WSS/Q: {wall_shear_stress(0, 0.5*H, num=10)/Q:.8e}")

    Q = G*L*H**3/(12.0*mu)*(1.0 - 0.630*H/L)
    print("Flow rate (approx H<<L): ", Q)

    # do a mesh plot of the velocity profiles
    x = np.linspace(0.0, L, 101)
    y = np.linspace(0.0, H, 101)
    xs, ys = np.meshgrid(x, y)
    zs = velocity(xs, ys, num=8)
    zs = zs/np.max(zs)
    print(np.min(zs))
    zs[zs<0] = 0

    fig, (ax, ax1) = plt.subplots(1, 2, constrained_layout=True, figsize=(10, max(1, 5.0/aspect_ratio)))
    cs = ax1.contourf(xs, ys, zs, 20, cmap=plt.cm.bone)
    cbar = plt.colorbar(cs)
    cbar.ax.set_ylabel("normalized velocity")

    wss = wall_shear_stress(0, y, num=8)
    ax.plot(wss/np.max(wss), y, "k-", lw=2)
    ax.set_xlabel("normalized WSS")
    ax.set_ylabel("height")

    fig.savefig(img + "_velocity" + ".pdf")


c1H = 0.5e-3
c1L = 3.2e-3
c2H = 0.8e-3
c2L = 1.2e-3

channel1_length_long = 16.4e-3 # big reservoir channel
channel1_length_short = 9.8e-3 # big reservoir channel
#channel_length = 9.8e-3 # small reservoir channel

# short channels
plot(aspect_ratio=c1L/c1H, channel_length=channel1_length_short, scale=c1L, img="channel1", mu=1e-3)
plot(aspect_ratio=c2L/c2H, channel_length=channel1_length_short, scale=c2L, img="channel2", mu=1e-3)

# long channels
plot(aspect_ratio=c1L/c1H, channel_length=channel1_length_long, scale=c1L, img="channel1", mu=1e-3)
plot(aspect_ratio=c2L/c2H, channel_length=channel1_length_long, scale=c2L, img="channel2", mu=1e-3)

# short channels
plot(aspect_ratio=c1L/c1H, channel_length=channel1_length_short, scale=c1L, img="channel1", mu=0.8e-3)
plot(aspect_ratio=c2L/c2H, channel_length=channel1_length_short, scale=c2L, img="channel2", mu=0.8e-3)

# long channels
plot(aspect_ratio=c1L/c1H, channel_length=channel1_length_long, scale=c1L, img="channel1", mu=0.8e-3)
plot(aspect_ratio=c2L/c2H, channel_length=channel1_length_long, scale=c2L, img="channel2", mu=0.8e-3)
