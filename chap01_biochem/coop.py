# -------------------------------------------------------------------

#  Python code for plotting the velocity of a cooperative enzyme.

#  For Chapter 1, Fig. 1.4 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Parameters
e0 = 1
k2 = 1
k4 = 2

s = np.linspace(0, 2, 100)

# Case 1
K1 = 1000
K2 = 0.001
V1 = (k2 * K2 + k4 * s) * e0 * s / (K1 * K2 + K2 * s + s ** 2)

# Case 2
K1 = 0.5
K2 = 2
V2 = (k2 * K2 + k4 * s) * e0 * s / (K1 * K2 + K2 * s + s ** 2)

# Case 3
K1 = 0.5
K2 = 100
V3 = (k2 * K2 + k4 * s) * e0 * s / (K1 * K2 + K2 * s + s ** 2)

# Plotting
plt.plot(s, V1, 'r', linewidth=2, label='K1=1000, K2=0.001')
plt.plot(s, V2, '--b', linewidth=2, label='K1=0.5, K2=2')
plt.plot(s, V3, '--g', linewidth=2, label='K1=0.5, K2=100')

plt.xlabel('s')
plt.ylabel('V')
plt.legend()
plt.show()
