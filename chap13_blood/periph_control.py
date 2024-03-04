
#   -------------------------------------------------------------------
# 
#    Make the plot for Fig. 13.10.
# 
#    For Chapter 13, Section 13.1.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

taulist = [3.8, 1.2]
alist = [0.36, 0.53]
klist = [3.15, 4.38]

omegabya = np.arange(0, np.pi/4, 0.01)

for j in range(2):
    tau = taulist[j]
    a = alist[j]
    kp1 = klist[j]
    omega = a * omegabya
    arg = omega * tau + kp1 * np.arctan(omegabya)
    alpha = -omega / np.tan(arg)
    mu = -omega / (np.cos(np.arctan(omegabya))**kp1 * np.sin(arg))
    ndx = np.where(alpha > 0)

    plt.figure(1)
    if j == 0:
        plt.plot(alpha[ndx], mu[ndx])
    if j == 1:
        plt.plot(alpha[ndx], mu[ndx], '--')

plt.axis([0, 4, -9, 0])
plt.text(0.5, -5, 'Unstable', fontsize=18)
plt.text(3, -2, 'Stable', fontsize=18)
plt.text(2.7, -5, 'Normal Human', fontsize=18)
plt.text(2.55, -8, 'CN Human', fontsize=18)

plt.plot([1.75, 1.75], [-9, 0], 'k', linewidth=0.5)
plt.plot([2.5, 2.5], [-9, 0], 'k', linewidth=0.5)
plt.xlabel(r'$\alpha$ (day$^{-1}$)')
plt.ylabel(r'$\mu$ (day$^{-1}$)')


