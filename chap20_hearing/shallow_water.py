
# -------------------------------------------------------------------
# 
#  Simulate shallow-water waves in the cochlea.
#  
#  For Chapter 20, Section 20.2.4 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
# 
# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

lam = 1.5
w_list = [800, 1500]

for w in w_list:
    L = 3.5
    l = 0.0035

    k0 = 1e7
    r0 = 3000

    alpha = np.sqrt(2 * w * w / (l * (k0 + 1j * w * r0)))
    ar = np.real(alpha)

    # Take the square root with positive imaginary part
    ai = -np.imag(alpha)

    x = np.linspace(0, L, 2000)
    eta = np.exp(3*lam*x/4 - 2*ai*np.exp(lam*x/2)/lam + 2*1j*ar*np.exp(lam*x/2)/lam)

    plt.figure(w_list.index(w) + 1)
    plt.plot(x, eta/np.max(np.abs(eta)), 'r')
    plt.plot(x, np.abs(eta)/np.max(np.abs(eta)), '--b')
    plt.plot(x, -np.abs(eta)/np.max(np.abs(eta)), '--b')
    xp = -2*np.log(4*ai/(3*lam))/lam
    plt.xlim([0, L])
    plt.xlabel('x (cm)')
    plt.title(f'$\omega = {w}$/s')
    plt.ylabel('normalized amplitude')

