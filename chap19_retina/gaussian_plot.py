
#  Code to plot excitatory and inhibitory surround function
 
#  For Figure 23, Chapter 19 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, s, g):
    return g * s * np.exp(-s**2 * x**2) / 2

# parameters
s1 = 3
s2 = 1
g1 = 3
g2 = 1
x = np.arange(-2, 2, 0.01)
G1 = gaussian(x, s1, g1)
G2 = gaussian(x, s2, g2)

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(x, G1, '--', label='g1(x)')
plt.plot(x, -G2, '--', label='g2(x)')
plt.plot(x, G1 - G2, label='f(x) = g1(x) + g2(x)')
plt.xlabel('x')
plt.legend()

