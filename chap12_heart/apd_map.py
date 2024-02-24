
#   -------------------------------------------------------------------
# 
#    APD map.
# 
#    For Chapter 12, Section 12.5.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Set parameters
A_max = 616
A = 750
mu = 170
DI_min = 100
BCL = 680
lim1 = BCL - DI_min
lim2 = 2 * BCL - DI_min

x1 = np.arange(0, lim1 + 1)
x2 = np.arange(lim1, lim2 + 1)
y1 = BCL - x1
y2 = 2 * BCL - x2

g1 = A_max - A * np.exp(-y1 / mu)
g2 = A_max - A * np.exp(-y2 / mu)

x = 400
nit = 120
g = np.zeros(2*nit)
# Find the alternans fixed point
for j in range(nit):
    y = BCL - x
    g[2*(j) - 1] = A_max - A * np.exp(-y / mu)
    g[2*(j)] = g[2*j - 1]
    x = g[2*(j) - 1]

plt.plot(x1, g1, x2, g2, np.concatenate([x1, x2]), np.concatenate([x1, x2]), '--', 
          435, 435, 'k*', 607, 607, 'k*', g[2*nit-8:2*nit-1], g[2*nit-7:2*nit], 'k', linewidth=2)
plt.axis([0, 1300, 0, 1000])
plt.xlabel(r'APD$_n$ (ms)', fontsize=12)
plt.ylabel(r'APD$_{n+1}$ (ms)', fontsize=12)

