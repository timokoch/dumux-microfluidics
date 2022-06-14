#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

angle = 18.0
theta = angle/180.0*np.pi
sinTheta = np.sin(theta)

rpm = 5.0
time = np.linspace(0, 60.0/rpm*3.0, 500)
t = time*2.0*np.pi*rpm/60.0
sint = np.sin(t)
cost = np.cos(t)
gamma = np.arcsin(-sinTheta*sint)
beta = np.arcsin(-sinTheta*cost/(np.sqrt(1.0 - sinTheta**2*sint**2)))

plt.plot(time, beta*180.0/np.pi)
plt.plot(time, gamma*180.0/np.pi)
plt.ylabel("Angle in degree")
plt.xlabel("Time in seconds")
plt.title(f"{rpm} rpm and {angle}° tilt")
plt.savefig("tilting_angles.pdf")
plt.show()

