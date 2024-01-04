#   -------------------------------------------------------------------
# 
#    Program to compute facilitation in a model of the residual bound
#    calcium hypothesis of synaptic plasticity.
# 
#    For Chapter 8, Fig. 8.9 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

tpcp = 200
k1 = 3.75e-3
k2 = 2.5e-3
k3 = 5e-4

km1 = 4e-4  # time units are ms
km2 = 1e-3
km3 = 0.1
K1 = km1 / k1
K2 = km2 / k2
K3 = km3 / k3

T = np.arange(1 / 0.186, 10000, 0.001)
a1 = np.exp(-km1 * (T + tpcp / K1))
a2 = np.exp(-km2 * (T + tpcp / K2))
a3 = np.exp(-km3 * (T + tpcp / K3))

Fmax = (1 / (1 - a1)) * (1 / (1 - a2)) * (1 / (1 - a3))

plt.semilogx(1000 / T, Fmax)  # the factor 1000 is to convert to Hertz = s^{-1}
plt.xlabel('Stimulus Frequency (Hz)')
plt.ylabel('Facilitation')
plt.axis([0.1, 100, 0, 7])
plt.show()
