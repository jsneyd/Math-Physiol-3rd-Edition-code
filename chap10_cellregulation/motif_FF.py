
#   -------------------------------------------------------------------
# 
#    The feedforward motif.
# 
#    For Chapter 10, Section 10.1.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
n = 2
K = 2
init = [0, 0]

# Define the ODE system
def rhs(t, Y):
    R, X = Y
    S = (t > 10) + (t > 20) + (t > 30)
    phi = K * X**n / (K + X**n)
    dRdt = S - phi * R
    dXdt = S - X
    return [dRdt, dXdt]

# Time points
tspan = np.linspace(0, 50, 5000)

# Integrate the ODEs
sol = solve_ivp(rhs, [tspan[0], tspan[-1]], init, t_eval=tspan)

# Extract results
T = sol.t
R = sol.y[0]
X = sol.y[1]
S = (T > 10).astype(int) + (T > 20).astype(int) + (T > 30).astype(int)

# Plot the results
plt.plot(T, R, label='R')
plt.plot(T, S, label='S')
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend()
plt.grid(True)
plt.show()
