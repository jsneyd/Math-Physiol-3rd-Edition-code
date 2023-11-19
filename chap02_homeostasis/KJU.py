# -------------------------------------------------------------------

#  Python code for solving the KJU model of a Na-transporting cell. The
#  nonlinear equations are solved using a built-in Python nonlinear solver, fsolve.

#  For Chapter 2, Fig. 2.21 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def KJUfun(x, Nm):
    Ks = 2.5
    Cs = 122.5
    Ns = Cs - Ks
    z = -2
    rhon = 1.5
    rhok = 1
    
    Ni, Ki, mu, v = x

    dum = np.exp(-v)
    Ci = Cs / dum

    out = [
        v * (Ni - Nm * dum) / (1 - dum) + 3 * rhon * Ni,
        v * (Ki - Ks * dum) / (1 - dum) - 2 * rhok * Ni,
        Ni + Ki - Ci + z / mu,
        Ni + Ki + Ci + 1 / mu - (Ns + Ks + Cs)
    ]

    return out


Nm = np.linspace(10, 210, 50)
x0 = [20, 100, 0.01, -3]
Ni = np.zeros_like(Nm)
Ki = np.zeros_like(Nm)
mu = np.zeros_like(Nm)
v = np.zeros_like(Nm)

for i in range(50):
    x = fsolve(lambda x: KJUfun(x, Nm[i]), x0)
    x0 = x
    Ni[i], Ki[i], mu[i], v[i] = x

plt.figure(1)
plt.plot(Nm, mu)
plt.xlabel('N_m')
plt.ylabel('scaled volume')

plt.figure(2)
plt.plot(Nm, v)
plt.xlabel('N_m')
plt.ylabel('scaled membrane potential')

plt.show()


