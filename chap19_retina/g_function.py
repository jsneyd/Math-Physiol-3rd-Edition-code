
#     -------------------------------------------------------------------
# 
#      Plot the function that controls light adaptation in the model of LA
#      in cones. Compare to data.
# 
#      For Chapter 19, Section 19.2.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

vstar = 35.7
s1 = 1.59 / vstar
s2 = 1130
vK = -13 / vstar
tauy = 0.07
k1 = 35.4
gam = 303
delta = 5
kappa = 0.1
eta = 52.5
tau1 = 0.012
taum = 0.016
tauz = 0.04

y = np.linspace(0.2, 1, 100)

phi = (y * np.exp(vK * (1 - y))) ** (1 / 3) * (delta + (gam - delta) * eta * (np.exp(-vK * (1 - y) / s1) - 1) / (s2 * k1 + eta * (np.exp(-vK * (1 - y) / s1) - 1)))
A = 4 + 84 / (1 + (y / 0.34) ** 4)

plt.plot(y, phi, 'r', label='Model')
plt.plot(y, A, '--b', label='Data')
plt.xlabel('y')
plt.ylabel('[cGMP]_dark s^-1')
plt.legend()
plt.show()
