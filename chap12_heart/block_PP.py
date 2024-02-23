
#   -------------------------------------------------------------------
# 
#    Phase portrait for blocking propagation failure.
# 
#    For Chapter 12, Section 12.4.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

R1 = 0.04
R2 = 1.0
a = 0.25

u = np.arange(0, 1.01, 0.01)

F = -u ** 4 / 4 + ((a + 1) * u ** 3) / 3 - a * u ** 2 / 2
V1 = np.sqrt(-F / R1)
F1 = F[-1]
V2 = np.sqrt((F1 - F) / R2)

plt.plot(u, V1, label='i=1')
plt.plot(u, V2, '--', label='i=2')
plt.xlabel('V', fontsize=20)
plt.ylabel('I', fontsize=20)
plt.legend()

