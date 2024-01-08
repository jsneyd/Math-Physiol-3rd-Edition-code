#   -------------------------------------------------------------------
# 
#    Code for the Hindmarsh-Rose model of electrical bursting.
# 
#    For Chapter 9, Section 9.1.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Global parameters
Iapp, r, s, x1 = 0.4, 0.001, 4, -(1 + np.sqrt(5)) / 2

# Simulation parameters
total = 1000
tstep = 0.01
Ilist = [0.4, 2, 4]
tspan = np.arange(0, total + tstep, tstep)
zstart = [0.2, 1.8, 2]

def deRHS(sol, t):
    x, y, z = sol
    zp = r * (s * (x - x1) - z)
    xp = y - x**3 + 3 * x**2 + Iapp - z
    yp = 1 - 5 * x**2 - y
    return [xp, yp, zp]

# Plotting
for j, Iapp in enumerate(Ilist):
    init = [0, 0, zstart[j]]
    sol = odeint(deRHS, init, tspan)

    plt.figure(2 * j + 1, figsize=(10, 6))
    plt.plot(tspan, sol[:, 0], label='x')
    plt.plot(tspan, sol[:, 2], '--', label='z')
    plt.title(f'I_app = {Iapp:.2f}', fontsize=18)
    plt.xlabel('t')
    plt.ylabel('x')
    plt.legend()
    plt.legend('boxoff')

    plt.figure(2 * j + 2, figsize=(10, 6))
    plt.plot(sol[:, 2], sol[:, 0])
    plt.title(f'I_app = {Iapp:.2f}', fontsize=18)
    plt.xlabel('z')
    plt.ylabel('x')

plt.show()
