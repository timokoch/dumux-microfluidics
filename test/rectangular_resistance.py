#!/usr/bin/env python3

import numpy as np


"""
Poiseuille flow circular tube (R)

dP/dx = 8µQ / (πR^4)
=> Q = - dP/dx πR^4 / (8µ)
     = - dP/dx K / µ
=> K = πR^4 / 8 = - Q/(dP/dx)*µ

u_max = 2*Q / πR^2

Poiseuille flow rectangular tube (H*L)
"""

# G = -dP/dx
G = 40/1e-2
mu = 1e-3 # Pa*s
H = 50e-6 # H in y-direction
L = 100e-6 # L in z-direction

def velocity(y, z, num):
    result = G/(2*mu)*y*(H - y)
    factor = -4.0*G*H**2/(mu*np.pi**3)
    for n in range(1, num+1):
        beta = (2.0*n - 1.0)*np.pi/H
        prefactor = 1.0/(2.0*n - 1.0)**3
        summand = factor*prefactor*(np.sinh(beta*z) + np.sinh(beta*(L - z)))/np.sinh(beta*L)
        result += summand
    return result


def flow_rate(num):
    result = G/(12*mu)*L*H**3
    factor = -16.0*G*H**4/(mu*np.pi**5)
    for n in range(1, num+1):
        beta = (2.0*n - 1.0)*np.pi/H
        prefactor = 1.0/(2.0*n - 1.0)**5
        summand = factor*prefactor * (np.cosh(beta*L) - 1.0)/np.sinh(beta*L)
        result += summand
    return result


def resistance(num):
    return flow_rate(num)/G*mu


for n in range(5, 15):
    print("Pressure gradient: ", G)
    print("Flow rate: ", flow_rate(n))
    print("Max velocity: ", velocity(y=0.5*H, z=0.5*L, num=n))

Q = G*L*H**3/(12.0*mu)*(1.0 - 0.630*H/L)
print("Flow rate (approx H<<L): ", Q)

R = np.sqrt(H*L/np.pi)
Q = G*np.pi*R**4/(8.0*mu)
print("Flow rate (equiv radius circular tube): ", Q)
print("Max velocity (equiv radius circular tube): ", 2.0*Q/(np.pi*R**2))
