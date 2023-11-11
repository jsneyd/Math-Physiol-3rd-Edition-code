# -------------------------------------------

# Python code for simulating the Goldbeter-Lefever model of glycolytic oscillations.

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
    f = u * (1 + w)**2
    Fu = par['nu'] - f
    Fw = f - par['eta'] * w

    return [Fu, Fw]


plt.close('all')
np.random.seed(42)  # For reproducibility

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['lines.linewidth'] = 2.0

# Parameters
par = {'nu': 190, 'eta': 120}

# Set up the integration
tspan = np.arange(0, 1.0005, 0.0005)  # time interval for the solution
IC = [37, 1]  # initial condition
sol = odeint(rhs, IC, tspan, args=(par,))

# Plot the solution as a function of time
plt.figure(1)
plt.plot(tspan, sol[:, 0], 'r', tspan, sol[:, 1], 'b--', linewidth=2)
plt.legend(['$\sigma_1$', '$\sigma_2$'], loc='upper right', fontsize=18)
plt.xlabel('time', fontsize=18)
plt.ylabel('concentration', fontsize=18)

# Plot the phase portrait with null clines
plt.figure(2)
w = np.arange(0.02, 20.1, 0.1)
u1 = par['nu'] / (1 + w)**2
u2 = par['eta'] * w / (1 + w)**2
plt.plot(sol[:, 0], sol[:, 1], 'r', u2, w, 'b--', u1, w, 'g--', linewidth=2)
plt.xlabel('$\sigma_1$', fontsize=20)
plt.ylabel('$\sigma_2$', fontsize=20)
plt.legend(['solution', '$\sigma_1$ nullcline', '$\sigma_2$ nullcline'], fontsize=18)
plt.axis([0, 40, 0, 20])

# Plot the Hopf bifurcation curve
eta = np.ones((201, 1)) * np.arange(0, 201)
nu = np.arange(0, 201).reshape(-1, 1) * np.ones((1, 201))
F = -eta**4 + eta**3 * nu - eta**3 - 3 * eta**2 * nu - 3 * eta * nu**2 - nu**3
plt.figure(3)
plt.contour(nu, eta, F, [0], linewidths=2)
plt.plot(par['nu'], par['eta'], '*')
plt.xlabel('$\\nu$')
plt.ylabel('$\eta$')
plt.text(50, 150, 'stable', fontsize=20)
plt.text(100, 60, 'unstable', fontsize=20)

plt.show()





