# -------------------------------------------------------------------
# 
#  Discrete release standing waves.
# 
#  For Chapter 7, Section 7.8 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
# 
# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

Abyk = 5
beta = np.arange(0.01, 5.1, 0.1)

astar = 1 - np.sqrt(1 / (1 + Abyk))
Af = Abyk * beta * np.sinh(beta)
l0 = np.cosh(beta)
lf = l0 + Af / 2
mu0 = l0 - np.sqrt(l0**2 - 1)
muf = lf - np.sqrt(lf**2 - 1)

C = Af / (Af + 2 * l0 - 2)

a = C * (muf - 1) / (muf - 1 / mu0)
b = C * (1 - 1 / mu0) / (muf - 1 / mu0)
bf = C - b * muf

plt.figure(1)
plt.plot(beta, a, 'r', label='a', linewidth=2)
plt.plot(beta, bf, 'b', label='bf', linewidth=2)
plt.xlabel(r'$\beta$', fontsize=16)
plt.ylabel(r'$c^*/c_e$', fontsize=16)
plt.text(2, 0.4, 'standing waves', fontsize=12)
plt.text(0.1, 0.1, 'traveling waves', fontsize=12)
plt.legend()


