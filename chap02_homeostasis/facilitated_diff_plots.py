# -------------------------------------------------------------------

#  Python code for computing fluxes and concentrations in a facilitated
#  diffusion system.

#  For Chapter 2, Fig. 2.3 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

# Parameters
rho = 10
s0 = 2
sL = 0.1
mu = 1 / ((1 + s0) * (1 + sL))

sig = np.arange(0, 20.1, 0.1)

u = sig / (1 + sig)

sig0 = s0 + rho * s0 / (1 + s0)
fac = (1 + mu * rho) * (s0 - sL)
xi = -(sig + rho * u - sig0) / fac
stot = sig + u
sigp = fac * (1 + sig) ** 2 / (sig ** 2 + rho + 2 * sig + 1)
up = sigp / (1 + sig) ** 2

# Plot 1
plt.figure(1)
plt.plot(xi, sig, 'r', xi, u, 'b', xi, stot, 'g', linewidth=2)
plt.axis([0, 1, 0, 3])
plt.xlabel('y', fontsize=20)
plt.ylabel('Oxygen Concentration', fontsize=20)
plt.text(0.1, 2.4, 'Total', fontsize=20)
plt.text(0.1, 0.45, 'Bound', fontsize=20)
plt.text(0.1, 1.3, 'Free', fontsize=20)
plt.show()

# Plot 2
plt.figure(2)
plt.plot(xi, sigp, xi, rho * up, '--', linewidth=2)
plt.text(0.3, 5, 'Bound Oxygen Flux', fontsize=20)
plt.text(0.3, 2.6, 'Free Oxygen Flux', fontsize=20)
plt.axis([0, 1, 0, 7])
plt.xlabel('y', fontsize=20)
plt.ylabel('Oxygen Flux', fontsize=20)
plt.show()

# Plot 3
rho = 5
gam = sig + rho * u
plt.figure(3)
plt.plot(gam, sig, sig, sig, '--', linewidth=2)
plt.axis([0, 10, 0, 10])
plt.xlabel('Oxygen Consumption', fontsize=20)
plt.ylabel('Critical external oxygen concentration', fontsize=20)
plt.text(6, 2.1, r'$\rho = 5$', fontsize=20)
plt.text(6, 5.8, r'$\rho = 0$', fontsize=20)
plt.show()

# Plot 4
r = rho
g = 14
xi = np.sqrt(gam / g)
xi0 = np.sqrt(sig / g)
s1 = g
s0 = (-r - s1 - 1 + np.sqrt(4 * r * s1 ** 2 + r ** 2 + 6 * r * s1 + s1 ** 2 + 2 * r + 2 * s1 + 1)) / (1 + s1) / 2
sig0 = np.arange(s0, 20.1, 0.1)
gam0 = sig0 + rho * sig0 / (1 + sig0)

xi1 = np.sqrt(np.abs(gam0 - s0 - r * s0 / (1 + s0)) / g)

plt.figure(4)
plt.plot(xi1, sig0, xi, sig, 'g', xi0, sig, 'r--', linewidth=2)
plt.xlabel('Radius', fontsize=20)
plt.ylabel('Free Oxygen', fontsize=20)
plt.legend([r'$\rho=5$', r'$\rho=5$, critical external concentration', r'$\rho=0$, critical external concentration'],
           loc='best')
plt.axis([0, 1, 0, 14])
plt.show()

