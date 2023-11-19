# -------------------------------------------------------------------

#  Python code for plotting the flux through a glucose transporter.

#  For Chapter 2, Fig. 2.7 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

def flux(se, si, k):
    return (se - si) / ((si + 1 + k) * (se + 1 + k) - k**2)

k = 0.5
se = np.arange(0, 10.1, 0.1)
si_values = [0, 1, 2]

flux_values = []

for si in si_values:
    flux_values.append(flux(se, si, k))

plt.figure(1)
plt.plot(se, flux_values[0], 'g', se, flux_values[1], 'r', se, flux_values[2], 'b', se, np.zeros_like(se), '--', linewidth=2)
plt.legend([r'$\sigma_i=0$', r'$\sigma_i=1$', r'$\sigma_i=2$'])
plt.xlabel('External glucose', fontsize=20)
plt.ylabel('Glucose flux', fontsize=20)
plt.show()
