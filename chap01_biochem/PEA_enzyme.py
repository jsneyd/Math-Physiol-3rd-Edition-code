#-------------------------------------------------------------------

# Python code for simulating the equilibrium approximation 
# of enzyme kinetics.

# For Chapter 1 of
# Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

# Written by James Keener and James Sneyd

#-------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def rhseq(y, t, a, b, e):
    s, z = y
    dsdt = z - s - b * s * (a + s - z)
    dzdt = e * (s - z)
    return [dsdt, dzdt]


# Part 1: Equilibrium approximation
a = 2.2
b = 1.7
e = 0.1

tspan = np.linspace(0, 20, 2000)
initial = [1, 1]
sol = odeint(rhseq, initial, tspan, args=(a, b, e))

s = sol[:, 0]
z = sol[:, 1]
x = (z - s) / a  # the original variable, x

# Plotting
plt.figure(1)
plt.plot(tspan, s, tspan, z, linewidth=2)
plt.legend(['$\sigma$', 'z'], fontsize=18)
plt.xlabel('$\\tau$', fontsize=20)
plt.axis([0, 2, 0, 1])

# Plot the slow manifold and the solution together
slow_s = np.linspace(0, 1, 100)
slow_z = slow_s + a * b * slow_s / (1 + b * slow_s)
slow_x = (slow_z - slow_s) / a

plt.figure(2)
plt.plot(s, z, slow_s, slow_z, '--', linewidth=2)
plt.xlabel('$\sigma$', fontsize=20)
plt.ylabel('z', fontsize=20)

plt.figure(3)
plt.plot(s, x, slow_s, slow_x, '--', linewidth=2)
plt.xlabel('$\sigma$', fontsize=20)
plt.ylabel('x', fontsize=20)

plt.show()


