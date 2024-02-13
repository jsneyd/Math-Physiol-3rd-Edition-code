
#   -------------------------------------------------------------------
# 
#    Code to solve the SNF circadian clock model of Kim and Forger.
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
    M, Pc, P = y
    A = (AT - P - Kd + ((AT - P - Kd)**2 + 4 * AT * Kd)**0.5) / 2
    Mp = alpha1 * A - bet1 * M
    Pcp = alpha2 * M - (bet2 + alpha3) * Pc
    Pp = alpha3 * Pc - bet3 * P
    return [Mp, Pcp, Pp]

# Parameters and initial conditions
Kd = 1e-3
alpha1 = 20
alpha2 = 1
alpha3 = 1
bet1 = 1
bet2 = 0
bet3 = 1
AT = 1
y0 = [0.5, 0.5, 0.5]  # Initial conditions

# Time points
t = np.linspace(-24, 30, 500)

# Solve the differential equations
sol = odeint(deRHS, y0, t)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, sol[:, 0], label='M')
plt.plot(t, sol[:, 1], label='Pc')
plt.plot(t, sol[:, 2], label='P')
plt.title('SNF Circadian Clock Model')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()

