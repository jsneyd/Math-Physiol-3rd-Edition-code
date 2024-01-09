
#   -------------------------------------------------------------------
# 
#    Code for the compound bursting model of Wierschem and Bertram, 2004.
#    Original code was provided by Richard Bertram and modified slightly.
# 
#    For Chapter 9, Section 9.1.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by Richard Bertram, James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

from numpy import exp,arange
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Global parameters
vn, sn, vm, sm, kd, kadp, gca, vca, gkca, vk, gkatp, gk, cm, f, alpha, kc, v1, tauc, v2, taun = [
    -16, 5.6, -20, 12, 0.3, 20, 1200, 25, 300, -75, 350, 3000, 5300, 0.001, 2.25e-6, 0.1, 10, 1.2e6, 185, 16
]

# Simulation parameters
total = 600000
tstep = 1

# Initial conditions
v0, n0, c0, atp0, adp0 = -61.21, 0.000311, 0.07, 1.99, 0.0502
u0 = [v0, n0, c0, atp0, adp0]

# Specify the output points
tspan = arange(0, total + tstep, tstep)

# ODE system for the compound bursting model
def deRHS(u, t):
    v, n, c, atp, adp = u

    # Activation variables
    ninf = 1 / (1 + exp((vn - v) / sn))
    minf = 1 / (1 + exp((vm - v) / sm))
    omega = 1 / (1 + (kd / c))

    # ATP handling
    phi = atp * (1 + kadp * adp)**2

    # Ionic currents
    ica = gca * minf * (v - vca)
    ikca = gkca * omega * (v - vk)
    ikatp = (gkatp / atp) * (v - vk)
    ik = gk * n * (v - vk)

    # Differential equations
    vp = (-ica - ik - ikatp - ikca) / cm
    np = (ninf - n) / taun
    cp = -f * (alpha * ica + kc * c)
    atpp = (v1 - phi) / tauc
    adpp = (phi - v2 * adp) / tauc

    return [vp, np, cp, atpp, adpp]

# Solve the ODE system
sol = odeint(deRHS, u0, tspan)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(tspan / 60000, sol[:, 0])
plt.ylabel('V')
plt.xlabel('t (s)')

plt.figure(figsize=(10, 6))
plt.plot(tspan / 1000, sol[:, 2])
plt.xlabel('t (s)')
plt.ylabel('Ca++')

plt.figure(figsize=(10, 6))
plt.plot(tspan / 60000, sol[:, 3])
plt.xlabel('t (s)')
plt.ylabel('ATP')

plt.show()
