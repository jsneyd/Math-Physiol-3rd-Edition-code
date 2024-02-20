
#   -------------------------------------------------------------------
# 
#    The repressilator (Elowitz and Leibler (2000) Nature, 403, 335).
# 
#    For Chapter 10, Section 10.1.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Parameters
alpha = 50
n = 3
m1_init = 0.7
m2_init = m1_init
m3_init = 0.8
init = [m1_init, m2_init, m3_init]
total_time = 50
tstep = 0.1

# Define the ODE system
def deRHS(sol, t):
    m1, m2, m3 = sol
    m1p = alpha / (1 + m2**n) - m1
    m2p = alpha / (1 + m3**n) - m2
    m3p = alpha / (1 + m1**n) - m3
    return [m1p, m2p, m3p]

# Time points
tspan = np.arange(0, total_time + tstep, tstep)

# Integrate the ODEs
sol = odeint(deRHS, init, tspan)

# Plot the results
plt.figure()
plt.plot(tspan, sol[:, 0], label='m1')
plt.plot(tspan, sol[:, 1], label='m2')
plt.plot(tspan, sol[:, 2], label='m3')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()
plt.grid(True)
plt.show()
