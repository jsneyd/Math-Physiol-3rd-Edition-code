# -------------------------------------------------------------------

#  Python code for plotting the heights of two water columns.

#  For Chapter 2, Fig. 2.15 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

h1 = np.arange(0, 2.01, 0.01)
p = h1**2 - h1
h2 = 2 - h1

plt.figure(1)
plt.plot(p, h1, p, h2)
plt.xlabel('p')
plt.legend([r'$\eta_1$', r'$\eta_2$'], loc='best')
plt.axis([0, 2, 0, 2])
plt.show()
