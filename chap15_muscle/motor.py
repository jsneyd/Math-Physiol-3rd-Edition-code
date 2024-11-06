
#   -------------------------------------------------------------------
# 
#    Velocity of a motor pulling a cargo.
# 
#    For Chapter 15, Section 15.10.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Define the function f(w0, v)
def getf(w0, v):
    wl = 10 * v
    return v - (wl ** 2) * (1 - np.exp(w0 - wl)) / \
           ((np.exp(w0) - 1) * (np.exp(-wl) - 1) + wl * (np.exp(w0 - wl) - 1))
           

# Generate values for w0 and v1
w0 = np.linspace(0, 10, 100)
D1 = 1
Dc_over_D1 = 0.1 / (1 + 0.1)
v1 = 2 * Dc_over_D1 * (np.exp(w0) - 1) / (np.exp(w0) + 1)

# Plot the contour and the dashed line
plt.plot(w0, v1, '--r')

n = 256
w0 = np.linspace(0, 10, n)
v = np.linspace(0, 0.4, n)
X, Y = np.meshgrid(w0, v)
plt.contour(X, Y, getf(X, Y),[0])

# Add labels and text
plt.xlabel(r'$\omega_0$')
plt.ylabel(r'$v\delta/D_1$')
plt.text(5, 0.16, 'hard spring')
plt.text(5, 0.33, 'soft spring')