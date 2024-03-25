
#   -------------------------------------------------------------------
# 
#    Mutual inhibition model of respiratory control.
# 
#    For Chapter 14, Section 14.7.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parameters for respiratory control
tau = 1.0
m = 0.5
mu = 5.0
k = 1.0
E1 = 0.5
E2 = 0.3
K = 0.2

# Define the right-hand side of the ODE system
def rhs(s, t):
    I1, I2, x, xdot = s

    arg1 = E1 - I2
    F1 = (2 * arg1 ** 2 / (K + arg1)) if arg1 > 0 else 0

    f = 2.5 * x ** 3 / (1 + x ** 3)

    arg2 = E2 - I1 + f
    F2 = (2 * arg2 ** 2 / (K + arg2)) if arg2 > 0 else 0

    FI1 = (F1 - I1) / tau
    FI2 = (F2 - I2) / tau
    Fx = xdot
    Fxd = (I1 - k * x - mu * xdot) / m

    return [FI1, FI2, Fx, Fxd]

# Initial conditions
init = [0.2, 0.2, 0.0, 0.0]

# Time vector
t = np.linspace(0, 100, 1000)

# Solve the ODE system
sol = odeint(rhs, init, t)

# Extract variables
I1, I2, x, xdot = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3]

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(t, I1, label='I1')
plt.plot(t, I2, '--', label='I2')
plt.xlabel('Time')
plt.ylabel('Activity')
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, E2 + 2.5 * x ** 3 / (1 + x ** 3))
plt.xlabel('Time')
plt.ylabel('E2 + f(x)')
plt.show()
