# -------------------------------------------------------------------

#  Python code for plotting the solutions of the pump-leak model.

#  For Chapter 2, Fig. 2.17 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# set parameter values
Ne = 437
Ke = 20
gamma = 0.11
zx = -1
RTbyF = 25.8

p = np.arange(0.01, 14.01, 0.01)  # pump rate values could go to 14.3

alpha = (Ne * np.exp(-3 * p) + Ke * np.exp(2 * p * gamma)) / (Ne + Ke)

mu = (1 + np.sqrt(-alpha * zx ** 2 + zx ** 2 + alpha)) / (2 * (1 - alpha))
y = (-zx + np.sqrt(zx ** 2 + 4 * alpha * mu ** 2)) / (2 * alpha * mu)

v = -RTbyF * np.log(y)
vna = RTbyF * (3 * p - np.log(y))
vk = RTbyF * (-2 * p * gamma - np.log(y))

plt.figure(1)
plt.plot(p, mu, 'r', linewidth=2)
plt.axis([0, 14, 0, 5])
plt.xlabel('Pump rate', fontsize=20)
plt.ylabel('Cell Volume', fontsize=20)

plt.figure(2)
plt.plot(p, v, 'r', p, vna, 'b', p, vk, 'g', linewidth=2)
plt.legend(['V', 'V_{Na}', 'V_K'], fontsize=18)
plt.xlabel('Pump rate P', fontsize=20)
plt.ylabel('Potential (mV)', fontsize=20)
plt.axis([0, 14, -100, 100])

# Now the modified model with P = rho u^3
u = y * np.exp(-3 * p)
rho = p / u ** 3

# plot everything as a function of rho
plt.figure(3)
plt.plot(rho, mu, 'r', linewidth=2)
plt.axis([0, 4, 0, 5])
plt.xlabel('Pump rate', fontsize=20)
plt.ylabel('Cell Volume', fontsize=20)

plt.figure(4)
plt.plot(rho, v, 'r', rho, vna, 'b', rho, vk, 'g', linewidth=2)
plt.legend(['V', 'V_{Na}', 'V_K'], fontsize=18)
plt.xlabel('Pump rate $\\rho$', fontsize=20)
plt.ylabel('Potential (mV)', fontsize=20)
plt.axis([0, 400, -80, 60])

plt.show()
