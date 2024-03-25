
#   -------------------------------------------------------------------
# 
#    Code to simulate post-inhibitory rebound coupled HH neurons.
# 
#    For Chapter 14, Section 14.7.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Global variables
gpir = 0.3
Vpir = 120
gL = 0.1
gsyn = 0.3
Vsyn = -80
phi = 3
thetasyn = -44
ksyn = 2
VL = -60

# Define the right-hand side of the ODE system
def rhs(s, t):
    V1, h1, V2, h2 = s

    minf1, hinf1, tauh1, Sinf1 = gates(V1)
    minf2, hinf2, tauh2, Sinf2 = gates(V2)

    V1p = -gpir * minf1**3 * h1 * (V1 - Vpir) - gL * (V1 - VL) - gsyn * Sinf2 * (V1 - Vsyn)
    h1p = phi * (hinf1 - h1) / tauh1
    V2p = -gpir * minf2**3 * h2 * (V2 - Vpir) - gL * (V2 - VL) - gsyn * Sinf1 * (V2 - Vsyn)
    h2p = phi * (hinf2 - h2) / tauh2

    return [V1p, h1p, V2p, h2p]

# Hodgkin-Huxley gate equations
def gates(V):
    minf = 1 / (1 + np.exp(-(V + 65) / 7.8))
    hinf = 1 / (1 + np.exp((V + 81) / 11))
    tauh = hinf * np.exp((V + 162.3) / 17.8)
    Sinf = 1 / (1 + np.exp(-(V - thetasyn) / ksyn))

    return minf, hinf, tauh, Sinf

# Initial conditions
init = [-50, 1, -40, 0]

# Time vector
t = np.arange(-50, 300, 0.1)

# Solve the ODE system
sol = odeint(rhs, init, t)

# Extract variables
V1, h1, V2, h2 = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3]

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(t, V1, label='V1')
plt.plot(t, V2, '--', label='V2')
plt.xlabel('Time')
plt.ylabel('Membrane Potential (mV)')
plt.legend()
plt.ylim(-80, -25)
plt.xlim(0, 300)
plt.show()
