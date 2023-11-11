#-------------------------------------------------------------------

# Python code for simulating the quasi-steady-state approximation (QSSA) 
# of enzyme kinetics.

# For Chapter 1 of
# Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

# Written by James Keener and James Sneyd

#-------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def rhsqss(y, t, a, k, e):
    s, x = y
    dsdt = -s + x * (s + a)
    dxdt = (s - x * (s + k)) / e
    return [dsdt, dxdt]



# Part 2: QSS approximation
a = 0.5
k = 1.5
e = 0.05

tspan = np.linspace(0, 20, 2000)
initial = [1, 0]
sol = odeint(rhsqss, initial, tspan, args=(a, k, e))

s = sol[:, 0]
x = sol[:, 1]

# Plotting
plt.figure(4)
plt.plot(tspan, s, tspan, x, linewidth=2)
plt.legend(['$\sigma$', 'x'], fontsize=18)
plt.xlabel('$\eta$', fontsize=20)
plt.axis([0, 4, 0, 1])

# Plot the slow manifold and the solution together
slow_s = np.linspace(0, 1.2, 100)
slow_x = slow_s / (slow_s + k)
plt.figure(5)
plt.plot(s, x, slow_s, slow_x, '--', linewidth=2)
plt.xlabel('$\sigma$', fontsize=20)
plt.ylabel('x', fontsize=20)
plt.axis([0, 1.2, 0, 0.5])

np.savetxt('temp3.dat', np.column_stack((tspan, s, x)), header='t s x')
np.savetxt('temp4.dat', np.column_stack((slow_s, slow_x)), header='s x')

plt.show()
