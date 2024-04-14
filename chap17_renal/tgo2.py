
#     -------------------------------------------------------------------
# 
#      Make \mu-\tau bifurcation curves in the Layton model of
#      tubuloglomerular oscillations.
# 
#      For Chapter 17, Section 17.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def Feval(tau):
    w = np.pi/(tau +k/(2*mu))
    return w/(2*np.sin(k*w/(2*mu)))-gamma
  

# Parameters
kbar = 0.5
cbar = np.exp(-kbar)
a = 10
kk = 0
k_values = np.arange(0.1, 2.1, 0.002)
U = np.zeros_like(k_values)
T = np.zeros_like(k_values)

# Iterate over k values
for idx, k in enumerate(k_values):
    c = np.exp(-k)
    F = 1 + np.tanh(a * (cbar - c))
    Fp = -a * (1/np.cosh(a * (cbar - c)))**2

    mu = k * F
    gamma = -k * Fp * c

    # Solve for tau
    tau = fsolve(lambda tau: Feval(tau), 0.2)[0]
    T[idx] = tau
    U[idx] = mu

# Plot mu-tau bifurcation curves
plt.figure()
plt.plot(U, T)
plt.xlabel('$\mu$')
plt.ylabel('$\tau$')
plt.axis([0, 1.2, 0, 2])
plt.text(0.45, 1, 'unstable', fontsize=18)
plt.text(0.45, 0.1, 'stable', fontsize=18)
plt.show()
