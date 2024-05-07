
#     -------------------------------------------------------------------
# 
#      Solutions of the Peskin lateral inhibition  model.
# 
#      For Chapter 19, Section 19.4.1 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

lam = 1
k = 20

t = np.arange(0, 3, 0.01)
R = 1/(lam + 1) - (k - 1)/(k - lam - 1) * np.exp(-k * t) + lam * k / ((k - lam - 1) * (lam + 1)) * np.exp(-(lam + 1) * t)

x = np.arange(-2, 2, 0.01)
E = (x >= 0)
sql = np.sqrt(lam + 1)
Ix = lam / (lam + 1) * ((x >= 0) * (1 - np.exp(-x * sql) / 2) + (x < 0) * np.exp(x * sql) / 2)

plt.figure()
plt.plot(t, R)
plt.xlabel('t')
plt.ylabel('R(t)')

plt.figure()
plt.plot(x, E, '--', label='stimulus, E')
plt.plot(x, E - Ix, label='response, R')
plt.xlabel('x')
plt.legend()



