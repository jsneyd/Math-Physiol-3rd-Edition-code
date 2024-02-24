
#   -------------------------------------------------------------------
# 
#    Plot reentry loop maps.
# 
#    For Chapter 12, Section 12.5.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Set parameters
L = 1
Tr = 3/4
t = np.arange(0, 5, 0.01)
dt = np.arange(Tr, 10, 0.01)
c = 1/2 + 2 * (dt - Tr) / (1 + (dt - Tr))
Lbyc = L / c

ndx = np.min(np.where(Lbyc <= Tr))

# Plot figure 1
plt.figure(1)
T = 1.75
plt.plot(dt[:ndx], Lbyc[:ndx], dt[ndx:], Lbyc[ndx:], ':', t, t, '--', t, Tr * np.ones(len(t)), ':', linewidth=3)
plt.plot([0, dt[ndx]], [T, T], ':', [dt[ndx], 5], [T, T], linewidth=3)
plt.axis([0, 2.5, 0, 2.5])
plt.xlabel(r'$\Delta T_n$')
plt.ylabel(r'$\Delta T_{n+1}$')
plt.text(0.2, 0.9, r'$T_r$', fontsize=20)
plt.text(0.2, 1.85, r'$T$', fontsize=20)
plt.text(0.8, 2.1, r'$L/c$', fontsize=20)

# Plot figure 2
plt.figure(2)
T = 1.25
plt.plot(dt[:ndx], Lbyc[:ndx], dt[ndx:], Lbyc[ndx:], ':', t, t, '--', t, Tr * np.ones(len(t)), ':', linewidth=3)
plt.plot([0, dt[ndx]], [T, T], ':', [dt[ndx], 5], [T, T], linewidth=3)
plt.axis([0, 2.5, 0, 2.5])
plt.xlabel(r'$\Delta T_n$')
plt.ylabel(r'$\Delta T_{n+1}$')
plt.text(0.2, 0.9, r'$T_r$', fontsize=20)
plt.text(0.2, 1.35, r'$T$', fontsize=20)
plt.text(0.8, 2.1, r'$L/c$', fontsize=20)

plt.show()
