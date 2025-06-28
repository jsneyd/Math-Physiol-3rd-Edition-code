
#   -------------------------------------------------------------------
# 
#    Code to solve the NNF circadian clock model of Kim and Forger.
# 
#    For Chapter 10, Section 10.3.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the differential equations
def deRHS(y, t):
    M, Pc, PT, AT, V = y
    A = (AT - PT - Kd + ((AT - PT - Kd)**2 + 4 * AT * Kd)**0.5) / 2
    Mp = alpha1 * A - bet1 * M
    Pcp = alpha2 * M - (bet2 + alpha3) * Pc
    PTp = alpha3 * Pc - bet3 * PT
    ATp = delta * (1 / (1 + V) - A)
    Vp = delta * (A - V)
    return [Mp, Pcp, PTp, ATp, Vp]

# Parameters and initial conditions

Kd = 1e-2
alpha1 = 20
alpha2 = 1
alpha3 = 1
bet1 = 1
bet2 = 0
bet3 = 1
AT = 1
delta = 0.5
y0 = [0.5, 0.5, 0.5, 1, 0.5]  # Initial conditions

# Time points
t = np.linspace(-100, 30, 500)

# Solve the differential equations
sol = odeint(deRHS, y0, t)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, sol[:, 0], label='M')
plt.plot(t, sol[:, 1], label='Pc')
plt.plot(t, sol[:, 2], label='PT')
plt.title('NNF Circadian Clock Model')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()
plt.grid(True)
plt.show()
