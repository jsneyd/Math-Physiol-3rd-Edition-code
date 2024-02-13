
#   -------------------------------------------------------------------
# 
#    The trp operon.
# 
#    For Chapter 10, Section 10.2.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Parameters
krat = 500
konrat = 5
krrat = 50
ktrat = 100
mbyk1 = 1
mbk2 = 2

# Time points
T = np.linspace(0, 50, 500)

# Calculate R and F
R = T**2 / (ktrat + T**2)
F = krat / (1 + konrat + krrat * R)

# Plot
plt.plot(T, F, label='F(T)')
plt.plot(T, mbyk1 * T, '--', label='$\mu/K=1$')
plt.plot(T, mbk2 * T, '--', label='$\mu/K=2$')

plt.xlabel('T')
plt.ylabel('Expression level')
plt.legend()
plt.grid(True)
plt.show()
