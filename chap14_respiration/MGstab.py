
#   -------------------------------------------------------------------
# 
#   Mackey-Glass stability plot.
# 
#    For Chapter 14, Section 14.6 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

y = np.arange(0, 3.01, 0.01)
n = 3

F = y**n / (1 + y**n)
Fp = n * y**(n - 1) / ((1 + y**n)**2)
g = F - y * Fp

plt.figure()
plt.plot(y, F, label='F(y)')
plt.plot(y, Fp, label='dF/dy')
plt.plot(y, g, label='g(y)')
plt.plot(y, np.zeros_like(y), '--', color='black')
plt.xlabel('y')
plt.legend()

