
#     -------------------------------------------------------------------
# 
#      Model of unregulated glomerular filtration.
# 
#      For Chapter 17, Section 17.1 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

# Define the function to find the root of
def rhs_root(Qe):
    Qi = (Pa - P1) / Ra
    dum = Qe / Qi + alpha * np.log((Qe / Qi - alpha) / (1 - alpha) + 0j) - 1 + KfL * pii / (alpha * Qi)
    return np.real(dum) 


# Use typical parameters to find values for the constants
P1 = 60
P2 = 18
Pa = 100
Pe = 0
Pd = 18
pii = 25  # mm Hg
Qi = 650
Qd = 125
Qe = Qi - Qd
alpha = pii / (P1 - P2)

# Calculate resistances
Ra = (Pa - P1) / Qi
Re = (P1 - Pe) / Qe
Rd = (P2 - Pd) / Qd
KfL = -(Qe / Qi + alpha * np.log((Qe / Qi - alpha) / (1 - alpha)) - 1) * (alpha * Qi) / pii


# note that I'm working backwards through the Pa values, and changing the initial condition
# for each iteration. And taking a LOT of iterations. This is because, if you don't get the initial
# condition close to the root, things break (due to other irrelevant roots). I couldn't get the 
# bisection method to work, either.

Pa_values = np.linspace(160,67,500)
# Iterate over Pa values and calculate Qe
Qe_values = []
init = 150
for Pa in Pa_values:
    Qe = optimize.fsolve(rhs_root, init)
    Qe_values.append(Qe[0])
    init = Qe[0]

# Plot results
Qi = (Pa_values - P1)/Ra
Qd = Qi - Qe_values
plt.plot(Pa_values, Qd, 'b', label='Qd')
plt.xlabel('Pa')
plt.ylabel('Qd')



