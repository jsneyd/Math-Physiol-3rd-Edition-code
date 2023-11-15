# -------------------------------------------------------------------

#  Python code for plotting Goldbeter-Koshland functions.

#  For Chapter 1, Fig. 1.6 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

w = np.linspace(0, 1, 100)

K1, K2 = 0.1, 0.05
vrat1 = (1 - w) * (K1 + w) / (w * (K2 + 1 - w))
plt.plot(vrat1, w, 'r', linewidth=2)

K1, K2 = 1.1, 1.2
vrat2 = (1 - w) * (K1 + w) / (w * (K2 + 1 - w))
plt.plot(vrat2, w, '--b', linewidth=2)

K1, K2 = 0.1, 1.2
vrat3 = (1 - w) * (K1 + w) / (w * (K2 + 1 - w))
plt.plot(vrat3, w, '--g', linewidth=2)

plt.xlabel('v1/v2')
plt.ylabel('w')
plt.xlim([0, 4])
plt.show()

