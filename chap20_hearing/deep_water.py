
# -------------------------------------------------------------------
# 
#  Simulate deep-water waves in the cochlea.
#  
#  For Chapter 20, Section 20.2.4 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
# 
# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

L = 3.5  # Units of cm
l = 3.5
N = 1000
x = np.linspace(0, L, N)
w_list = [800, 1500]
igorout = []

for w in w_list:
    sigma = l / L
    k0 = 1e7
    r0 = 3000
    F0 = 1000
    lam = 1.5
    rho = 1
    k = k0 * np.exp(-lam * x)
    r = r0 * np.exp(-lam * x)
    Z = (r + k / (1j * w))
    Z0 = (r0 + k0 / (1j * w))
    Y = 2 / Z
    eta0 = F0 / Z0

    intY = 2j * w * (np.exp(lam * x) - 1) / (lam * (1j * w * r0 + k0))
    eta = -1j * rho * eta0 * Y * np.exp(-w * rho * intY)

    plt.figure(w_list.index(w) + 1)
    plt.plot(x, eta / np.max(np.abs(eta)), 'r')
    plt.plot(x, np.abs(eta) / np.max(np.abs(eta)), '--b', linewidth=2)
    plt.plot(x, -np.abs(eta) / np.max(np.abs(eta)), '--b', linewidth=2)
    plt.xlim([0, L])
    plt.title(f'$\omega = {w}$/s')
    plt.ylabel('normalized amplitude')
    plt.xlabel('x (cm)')

    igorout.append([x, np.real(eta / np.max(np.abs(eta))), np.abs(eta / np.max(eta)), -np.abs(eta / np.max(eta))])

# writematrix(igorout,'deep.dat')  # for external plotting
plt.show()
