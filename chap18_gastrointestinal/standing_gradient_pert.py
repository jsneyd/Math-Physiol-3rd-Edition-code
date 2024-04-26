
#     -------------------------------------------------------------------
# 
#  Plot the singular perturbation solution of the standing gradient osmotic
#  solution.
# 
#      For Chapter 18, Section 18.3.1 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# parameters
D = 1000
r = 0.05
c0 = 0.3
P = 0.2
L = 100

alp = 0.1

N0 = 0.3

eps = D * r / (L**2 * c0 * P)
n0 = N0 / (c0**2 * P)
ud = 1/2 + np.sqrt(1 + 4 * n0) / 2
wa = 2 * alp * (ud - 1)

# the outer solution
wo = wa * ud * np.linspace(0, 1, 1000)
uo = wa * ud / wo
y = (wa - wo - wa * ud * np.log((wo - wa * ud) / (wa - wa * ud))) / 2 + alp

# Composite solution
ndx = np.where((y >= alp) & (y <= 1))
w1 = wo[ndx][-1]
u1 = uo[ndx][-1]
tmp = (1 - u1) * np.exp(w1 * (y[ndx] - 1) / eps)
Y = np.concatenate(([0, alp], y[ndx]))
Uo = np.concatenate(([ud, ud], uo[ndx] + tmp))
plt.plot(y[ndx], uo[ndx], '--', Y, Uo)
plt.legend(['outer solution', 'composite solution'], loc='best')
plt.xlabel('y')
plt.ylabel('u')
plt.show()
