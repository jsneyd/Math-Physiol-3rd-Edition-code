#   -------------------------------------------------------------------
# 
#    Code for the fast subsystem of the Hindmarsh-Rose model of electrical bursting.
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

# Function to calculate the right-hand side of the ODEs
def deRHS(s, t):
    x, y = s
    xp = y - x**3 + 3 * x**2 + Iapp
    yp = 1 - 5 * x**2 - y
    return [xp, yp]

# change Iapp to get different phase portraits
Iapp = 0

# Phase plane
x = np.arange(-2, 2.01, 0.01)
y1 = x**3 - 3 * x**2 - Iapp
y2 = 1 - 5 * x**2
x1, x2 = -1/2 - np.sqrt(5)/2, -1/2 + np.sqrt(5)/2
ys1, ys2 = 1 - 5 * x1**2, 1 - 5 * x2**2

# Initial data
init = [0.62, -0.86]

# Simulation parameters
total = 100
tstep = 0.01
tspan = np.arange(0, total + tstep, tstep)

# Solve the ODE
sol = odeint(deRHS, init, tspan)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(sol[:, 0], sol[:, 1], label='solution trajectory')
plt.plot(x, y1, '--', label='dx/dt=0')
plt.plot(x, y2, '--', label='dy/dt=0')
plt.axis([-2, 2, -8, 2])
plt.legend()
plt.xlabel('x')
plt.ylabel('y')

if Iapp == 0:
    plt.plot(x1, ys1, '*', markersize=12)
    plt.plot(x2, ys2, '*', markersize=12)

plt.show()
