# -------------------------------------------

# Python code for simulating the Sel'kov model of glycolytic oscillations.

# For Chapter 1 of
# Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

# Written by James Keener and James Sneyd.

# -------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def rhs(y, t, par):
    u = y[0]
    w = y[1]
    g = par['g']
    f = u * w**g / (w**g * u + w**g + 1)
    Fu = par['nu'] - f
    Fw = par['alph'] * f - par['eta'] * w

    return [Fu, Fw]


plt.close('all')
np.random.seed(42)  # For reproducibility

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['lines.linewidth'] = 2.0

par = {'nu': 0.0285, 'eta': 0.08, 'alph': 1.0, 'g': 2}

# Set up the integration
tspan = np.arange(0, 800.1, 0.1)  # time interval for the solution
IC = [0.4, 0.4]  # initial condition
sol = odeint(rhs, IC, tspan, args=(par,))

# Plot the solution as a function of time
plt.figure(1)
plt.plot(tspan, sol[:, 0], 'r', tspan, sol[:, 1], 'b--', linewidth=2)
plt.legend(['$\sigma_1$', '$\sigma_2$'], loc='upper right', fontsize=18)
plt.xlabel('time', fontsize=18)
plt.ylabel('concentration', fontsize=18)

# Plot the phase portrait
plt.figure(2)

# Add null clines
s2 = np.arange(0.02, 1.4, 0.01)
s1a = par['nu'] / (1 - par['nu']) * (1 + s2**par['g']) / s2**par['g']
s1b = (1 + s2**par['g']) / (s2**(par['g'] - 1) * (par['alph'] / par['eta'] - s2))
plt.plot(sol[:, 0], sol[:, 1], 'r', s1a, s2, 'b--', s1b, s2, 'g--', linewidth=2)
plt.legend(['solution', '$\sigma_1$ nullcline', '$\sigma_2$ nullcline'], fontsize=20)
plt.xlabel('$\sigma_1$', fontsize=20)
plt.ylabel('$\sigma_2$', fontsize=20)
plt.axis([0, 1.4, 0, 1])

eta = np.ones((51, 1)) * np.arange(0, 1.01, 0.01)
nu = np.arange(0, 0.51, 0.01).reshape(-1, 1) * np.ones((1, 101))

F = 2 * eta**3 * nu + nu**4 - eta**3 + eta * nu**2 - 2 * nu**3 + nu**2
plt.figure(3)
plt.contour(nu, eta, F, [0], linewidths=2)
plt.plot(par['nu'], par['eta'], '*')
plt.xlabel('$\\nu $')
plt.ylabel('$\eta$')
plt.text(0.1, 0.6, 'unstable', fontsize=20)
plt.text(0.3, 0.2, 'stable', fontsize=20)

plt.show()

