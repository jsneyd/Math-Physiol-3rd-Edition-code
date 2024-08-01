
#   -------------------------------------------------------------------
# 
#   Plot solutions of the lac operon model 
# 
#    For Chapter 10, Section 10.2.5 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Initial conditions and time span
init = [0.01, 0.01, 0.01, 0.0, 0.0]
tspan = np.linspace(0, 2000, 1000)

# Parameters for ODE solver
tol = 1e-9

# Function to calculate Le
def getLe(t):
    if t < 1000:
        return 0.04 * t / (50 + t)
    else:
        return 0.04 * 1000 / (50 + 1000) * 1 / (1 + (t - 1000) / 500)

# Function for the system of ODEs
def rhs(t, x):
    M, B, P, A, L = x
    
    aA, bA, gamA = 1.76e4, 2.15e4, 0.52
    aB, gamB = 1.66e-2, 2.26e-2
    aP, gamP = 10.0, 0.65
    aM, gamM = 9.97e-4, 0.41
    aL, gamL = 2880.0, 2.26e-2
    Kaa, K, KA, KL, KLe, KL1 = 2.52e4, 6000.0, 1.95, 9.7e-7, 0.26, 1.81
    
    Le = getLe(t)
    
    dM = aM * (1.0 + Kaa * A**2) / (K + Kaa * A**2) - gamM * M
    dB = aB * M - gamB * B
    dP = aP * M - gamP * P
    dA = aA * B * L / (KL + L) - bA * B * A / (KA + A) - gamA * A
    dL = aL * P * Le / (KLe + Le) - aA * B * L / (KL + L) - gamL * L
    
    return [dM, dB, dP, dA, dL]

# Solve the system of ODEs
sol = solve_ivp(rhs, [tspan[0], tspan[-1]], init, method='LSODA', t_eval=tspan, atol=tol, rtol=tol, max_step=0.01)

# Calculate Le for each time point
Le = np.array([getLe(t) for t in sol.t])

# Plot the results
plt.figure(1)
plt.plot(sol.t, sol.y[3], label='A', linewidth=2)
plt.ylabel('A')
plt.xlabel('Time')

plt.twinx()
plt.plot(sol.t, Le, label='Le', linewidth=2, color='orange')
plt.ylabel('Le')
plt.legend()

