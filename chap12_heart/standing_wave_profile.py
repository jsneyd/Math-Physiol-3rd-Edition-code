
#   -------------------------------------------------------------------
# 
#    Calculate the potentials in a passive cardiac cable.
# 
#    For Chapter 12, Section 12.4.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

L = 0.012
qi = 5.47e-3
qe = 0.5 * qi
Q = qi + qe
sQ = np.sqrt(Q)
lambda_ = sQ / L
lambda_g = 0.09
mu = np.exp(-L / lambda_g)
m = mu
ratc = qi / Q
E = np.exp(sQ)
LAM = 2 * (mu - E) * (mu - 1 / E) / (mu * (E - 1 / E))

A = -E / (E ** 2 * LAM + 2 * E ** 2 - 4 * E * m - LAM + 2)
B = A
C = -(2 * E ** 2 * LAM * m ** 2 + 2 * E ** 2 * m ** 2 * ratc - 2 * E ** 2 * LAM * m + 2 * E ** 2 * m ** 2 - 4 * E * m ** 3
      - 4 * E * m ** 2 * ratc - 4 * E ** 2 * m - 2 * E ** 2 * ratc + 8 * E * m ** 2 - 2 * LAM * m ** 2
      + 2 * m ** 2 * ratc + 2 * E ** 2 - 4 * E * m + 4 * E * ratc + 2 * LAM * m + 2 * m ** 2 - 4 * m - 2 * ratc + 2) \
     / ((E ** 2 * LAM + 2 * E ** 2 - 4 * E * m - LAM + 2) * (m ** 2 - 2 * m + 1))

x = np.arange(0, 1.05, 0.05) * L

phi = (mu - 1 / E) * np.exp(lambda_ * x) + (mu - E) * np.exp(-lambda_ * x)
Pk = 2 * qe / Q * (mu - 1 / E) * (mu - E) / (mu - 1)
phi_i = phi * qi / Q + Pk
phi_e = -phi * qe / Q + Pk

phin = (1 / mu - 1 / E) * np.exp(lambda_ * x) + (1 / mu - E) * np.exp(-lambda_ * x)
Pkn = 2 * qe / Q * (1 / mu - 1 / E) * (1 / mu - E) / (1 / mu - 1)
phi_in = (phin * qi / Q + Pkn)
phi_en = -phin * qe / Q + Pkn

Vi = []
Ve = []
y = []

for j in range(-25, 0):
    Vi = np.append(Vi, B * mu ** (-j) * phi_in + 1 + C)
    Ve = np.append(Ve, B * mu ** (-j) * phi_en + C)
    y = np.append(y, j * L + x)

for j in range(26):
    Vi = np.append(Vi, A * mu ** j * phi_i)
    Ve = np.append(Ve, A * mu ** j * phi_e)
    y = np.append(y, j * L + x)

plt.figure(1)
plt.plot(y, Ve, label='Ve', linewidth=2)
plt.plot(y, Vi, label='Vi', linewidth=2)
plt.plot(y, Vi - Ve, label='Vi - Ve', linewidth=2)
plt.xlabel('Length (cm)', fontsize=18)
plt.ylabel('Potential', fontsize=18)
plt.axis([-0.3, 0.3, -0.2, 1])
plt.legend()

