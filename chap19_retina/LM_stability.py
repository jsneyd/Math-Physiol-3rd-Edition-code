
#  --------------------------------
#  Code to plot the stability curve for  the Longtin-Milton model of the pupil
#  light reflex.  

#  Used to generate the image in Fig. 26 of Chapter 19 
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd
#  --------------------------------

import numpy as np
import matplotlib.pyplot as plt

eta = np.arange(np.pi / 2 + 0.01, np.pi, 0.01)
G = -1 / np.cos(eta)
delay = -eta / np.tan(eta)

plt.figure(figsize=(8, 6))
plt.plot(delay, G)
plt.axis([0, 10, 0, 10])
plt.xlabel('Delay')
plt.ylabel('Gain')
plt.text(1, 1, 'Stable', fontsize=18)
plt.text(4, 4, 'Unstable', fontsize=18)
plt.show()
