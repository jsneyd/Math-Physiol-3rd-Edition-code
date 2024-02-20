
#   -------------------------------------------------------------------
# 
#    Calculate sawtooth potentials in a passive cardiac cable.
# 
#    For Chapter 12, Section 12.4.1 of
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

E = np.exp(sQ)

x = np.arange(0, 1.05, 0.05) * L

phi = (mu - 1 / E) * np.exp(lambda_ * x) + (mu - E) * np.exp(-lambda_ * x)

Pk = 2 * qe / Q * (mu - 1 / E) * (mu - E) / (mu - 1)
phi_i = phi * qi / Q + Pk
phi_e = -phi * qe / Q + Pk
Vi = []
Ve = []
y = []

sc = 1.e2
for j in range(16):
    Vi = np.append(Vi, sc * mu ** j * phi_i)
    Ve = np.append(Ve, sc * mu ** j * phi_e)
    y = np.append(y, j * L + x)

plt.figure(figsize=(10, 6))
plt.plot(y, Ve, label='Ve', linewidth=2)
plt.plot(y, Vi, label='Vi', linewidth=2)
plt.plot(y, Vi - Ve, label='Vi - Ve', linewidth=2)
plt.xlabel('Length (cm)', fontsize=18)
plt.ylabel('Potential', fontsize=18)
plt.text(0.05, 0, 'Ve', fontsize=18)
plt.text(0.05, -18, 'V', fontsize=18)
plt.text(0.05, -10, 'Vi', fontsize=18)
plt.legend()

